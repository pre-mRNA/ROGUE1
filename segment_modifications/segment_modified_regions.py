import time
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from multiprocessing import Pool

def parse_bed_file(filepath):
    
    total_sites = 0
    parsed_sites = 0
    start_time = time.time()

    data = defaultdict(lambda: defaultdict(list))

    min_coverage = 10

    with open(filepath, 'r') as file:
        for line in file:
            total_sites += 1
            fields = line.strip().split()
            if len(fields) < 11:
                continue

            # parse modkit pileup outoput 
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            coverage = int(fields[4])
            strand = fields[5]
            stoichiometry = float(fields[10])

            # only consider cites with sufficient coverage 
            if coverage >= min_coverage:
                parsed_sites += 1
                data[chrom][strand].append((start, end, stoichiometry))

    elapsed_time = time.time() - start_time

    print(f"Total sites in input: {total_sites}")
    print(f"Total sites parsed with coverage at least {min_coverage}: {parsed_sites}")
    print(f"Processing time: {elapsed_time:.2f} seconds")

    return data

def process_chromosome_strand(data_tuple):
    chrom, strand, sites = data_tuple
    results = defaultdict(list)
    window_size = 100
    step_size = 10
    m6a_threshold = 10.0  # m6A site stoichiometry threshold 
    required_m6a_sites = 3
    required_tested_sites = 10
    
    m6A_count = 0  
    window_types = {'hypermethylated': 0, 'non-m6A': 0, 'intermediate': 0}

    # use intervaltree to segment sites 
    tree = IntervalTree(Interval(start, start + 1, (start, end, stoichiometry)) for start, end, stoichiometry in sites)
    sites.sort()
    last_position = max(sites, key=lambda x: x[0])[0] if sites else 0

    # using a sliding window approach 
    for i in range(sites[0][0], last_position, step_size):
        window_end = i + window_size
        window_sites = [site.data for site in tree[i:window_end]]
        if len(window_sites) < required_tested_sites:
            continue  # skip windows with insufficient tested sites 

        m6a_count = sum(1 for s in window_sites if s[2] > m6a_threshold)
        m6A_count += m6a_count  

        unique_sites = set(s[0] for s in window_sites)  # unique site identification 
        if m6a_count >= required_m6a_sites:
            cluster_type = 'hypermethylated'
            window_types['hypermethylated'] += 1
            color = '255,0,0' if strand == '+' else '0,0,255'

        elif m6a_count == 0:
            cluster_type = 'non-m6A'
            window_types['non-m6A'] += 1
            color = '128,128,128' if strand == '+' else '192,192,192'
 
        else:
            cluster_type = 'intermediate'
            window_types['intermediate'] += 1
            color = '255,255,0' if strand == '+' else '0,255,255'


        results[chrom].append({
            'start': i,
            'end': window_end,
            'type': cluster_type,
            'score': m6a_count,
            'strand': strand,
            'color': color,
            'unique_sites': unique_sites
        })

    return results, m6A_count, window_types  

def identify_clusters(data):
    pool = Pool(processes=104)  

    data_tuples = [(chrom, strand, sites) for chrom, strands in data.items() for strand, sites in strands.items()]
    results = pool.map(process_chromosome_strand, data_tuples)

    cluster_results = defaultdict(lambda: defaultdict(list))
    total_m6A_count = 0
    total_window_types = {'hypermethylated': 0, 'non-m6A': 0, 'intermediate': 0}

    # unpack results 
    for result, m6A_count, window_types in results:
        for chrom in result:
            for cluster in result[chrom]:
                cluster_results[chrom][cluster['strand']].append(cluster)
        total_m6A_count += m6A_count
        for window_type, count in window_types.items():
            total_window_types[window_type] += count

    pool.close()
    pool.join()

    # print summary 
    print(f"Total number of m6A sites found: {total_m6A_count}")
    for window_type, count in total_window_types.items():
        print(f"Total {window_type} windows: {count}")

    return cluster_results, total_m6A_count

# function to subtract m6A ranges 
def handle_subtraction(chrom, strands, hypermethylated, non_m6A):
    def subtract_intervals(base_intervals, subtracting_intervals):
        subtract_tree = IntervalTree(Interval(i['start'], i['end']) for i in subtracting_intervals)
        result = []
        for base in base_intervals:
            base_start, base_end = base['start'], base['end']
            overlaps = subtract_tree[base_start:base_end+1]
            points = {base_start, base_end + 1}
            for ov in overlaps:
                points.add(ov.begin)
                points.add(ov.end)
            sorted_points = sorted(points)
            last_point = base_start
            for point in sorted_points:
                if point > last_point and not any(ov.begin <= last_point < ov.end for ov in overlaps):
                    result.append({
                        'start': last_point,
                        'end': point - 1,
                        'type': base['type'],
                        'score': base['score'],
                        'strand': base['strand'],
                        'color': base.get('color', '0,0,0') 
                    })
                last_point = point
        return result

    result = []
    for strand, intervals in strands.items():
        hyper_intervals = [i for i in hypermethylated if i['strand'] == strand]
        non_m6A_intervals = [i for i in non_m6A if i['strand'] == strand]
        subtracted_intervals = subtract_intervals(non_m6A_intervals, hyper_intervals)
        for interval in subtracted_intervals:
            interval['strand'] = strand 
        result.extend(subtracted_intervals + [i for i in intervals if i['type'] != 'non-m6A' and 'strand' in i])
    return chrom, result
def merge_clusters(clusters, total_m6A_count):
    priority = {'hypermethylated': 3, 'intermediate': 2, 'non-m6A': 1}
    
    tasks = [(chrom, strands, [i for i in strands.values() for i in i if i['type'] == 'hypermethylated'], 
              [i for i in strands.values() for i in i if i['type'] == 'non-m6A']) for chrom, strands in clusters.items()]

    merged_clusters = defaultdict(lambda: defaultdict(list))
    with Pool() as pool:
        results = pool.starmap(handle_subtraction, tasks)
        for chrom, strand_results in results:
            for res in strand_results:
                merged_clusters[chrom][res['strand']].append(res)

    for chrom, strands in merged_clusters.items():
        for strand, intervals in strands.items():
            sorted_intervals = sorted(intervals, key=lambda x: (x['start'], -priority[x['type']]))
            merged = []
            for interval in sorted_intervals:
                if not merged or interval['start'] > merged[-1]['end'] or priority[interval['type']] > priority[merged[-1]['type']]:
                    merged.append(interval)
                else:
                    merged[-1]['end'] = max(merged[-1]['end'], interval['end'])
                    merged[-1]['score'] += interval['score']  # merge score for overlapping ranges 
            merged_clusters[chrom][strand] = merged

    # validate counts 
    merged_m6A_count = sum(
        interval['score'] for strand_intervals in merged_clusters.values() 
        for intervals in strand_intervals.values() for interval in intervals
    )
    if merged_m6A_count != total_m6A_count:
        print(f"Error: Merged m6A site count ({merged_m6A_count}) does not match initial count ({total_m6A_count}).")

    return merged_clusters


def output_bed_format(clusters, output_file_path):
    summary = {
        'hypermethylated': {'count': 0, 'total_width': 0},
        'intermediate': {'count': 0, 'total_width': 0},
        'non-m6A': {'count': 0, 'total_width': 0},
    }

    with open(output_file_path, 'w') as file:
        
        # file.write("Cluster details in BED format:\n")
        for chrom, strands in clusters.items():
            for strand, cluster_list in strands.items():
                for cluster in cluster_list:
                    region_type = cluster['type']
                    region_width = cluster['end'] - cluster['start']
                    summary[region_type]['count'] += 1
                    summary[region_type]['total_width'] += region_width
                    file.write(f"{chrom}\t{cluster['start']}\t{cluster['end']}\t{cluster['type']}\t{cluster['score']}\t{strand}\t0\t0\t{cluster['color']}\n")

    print("Summary of clustering results:")
    for region_type, stats in summary.items():
        if stats['count'] > 0:
            average_width = stats['total_width'] / stats['count']
        else:
            average_width = 0
        print(f"{region_type.capitalize()} regions: {stats['count']}, Total width: {stats['total_width']}, Average width: {average_width:.2f} nt")

# for testing 
def output_clusters_before_merging(clusters, output_file_path):
    with open(output_file_path, 'w') as file:
        for chrom, strands in clusters.items():
            for strand, cluster_list in strands.items():
                for cluster in cluster_list:
                    file.write(
                        f"{chrom}\t{cluster['start']}\t{cluster['end']}\t{cluster['type']}\t{cluster['score']}\t{strand}\t0\t0\t{cluster['color']}\n"
                    )



# in and out files 
bed_file = "/g/data/lf10/as7425/2023_splicing_termination/analysis/2024-06-15_intron-m6A-CNN//intron_nascent_modkit_0.995_sites.bed"
m6A_regions = "/g/data/lf10/as7425/2023_splicing_termination/analysis/2024-06-15_intron-m6A-CNN/hypermethylated_regions_v2.bed"
m6A_regions_before_merging = "/g/data/lf10/as7425/2023_splicing_termination/analysis/2024-06-15_intron-m6A-CNN/hypermethylated_regions_before_merging.bed"

# read in modification pileup 
data = parse_bed_file(bed_file)

# create clusters 
# each modification site is uniquely assigned to a 0-width cluster
# clusters are spaced at {interval_size} (currently 10 nt)
clusters, total_m6A_count = identify_clusters(data)

# print the intermediate clusters 
output_clusters_before_merging(clusters, m6A_regions_before_merging)

merged_clusters = merge_clusters(clusters, total_m6A_count)

output_bed_format(merged_clusters, m6A_regions)