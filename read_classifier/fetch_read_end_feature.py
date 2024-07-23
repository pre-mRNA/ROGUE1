import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import os
import tempfile

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

def get_end_feature(temp_bed_files, genome_file, output_dir, dog_bed_file, exon_bed_file, intron_bed_file, gene_id_map, min_overlap=5):
    
    logging.info(f"Processing {len(temp_bed_files)} temp bed files")
    
    def process_bed_chunk(bed_chunk):
        base_name = os.path.splitext(os.path.basename(bed_chunk))[0]
        with tempfile.NamedTemporaryFile(mode='w+t', dir=output_dir, suffix=f'_{base_name}_3prime_15nt.bed', delete=False) as output_file:
            output_path = output_file.name
        
        # extract 3' 15 nt and filter invalid reads
        process_cmd = f"""
        awk 'BEGIN {{OFS="\t"}}
        {{
            if (!($4 in chrom)) {{
                chrom[$4] = $1;
                strand[$4] = $6;
                total_width[$4] = 0;
                intervals[$4] = "";
            }}
            if (chrom[$4] != $1 || strand[$4] != $6) {{
                print $4 > "{output_dir}/{base_name}_invalid_reads.txt";
                next;
            }}
            total_width[$4] += $3 - $2;
            intervals[$4] = intervals[$4] " " $0;
        }}
        END {{
            for (read in total_width) {{
                if (total_width[read] < 15) {{
                    print read > "{output_dir}/{base_name}_invalid_reads.txt";
                    continue;
                }}
                split(intervals[read], arr);
                if (strand[read] == "+") {{
                    remaining = 15;
                    for (i = NF; i > 0; i -= 6) {{
                        start = arr[i-4];
                        end = arr[i-3];
                        width = end - start;
                        if (width >= remaining) {{
                            print arr[i-5], end - remaining, end, read, arr[i-1], arr[i];
                            break;
                        }} else {{
                            print arr[i-5], start, end, read, arr[i-1], arr[i];
                            remaining -= width;
                        }}
                    }}
                }} else {{
                    remaining = 15;
                    for (i = 1; i <= NF; i += 6) {{
                        start = arr[i+1];
                        end = arr[i+2];
                        width = end - start;
                        if (width >= remaining) {{
                            print arr[i], start, start + remaining, read, arr[i+4], arr[i+5];
                            break;
                        }} else {{
                            print arr[i], start, end, read, arr[i+4], arr[i+5];
                            remaining -= width;
                        }}
                    }}
                }}
            }}
        }}' {bed_chunk} | sort -k1,1 -k2,2n > {output_path}
        """
        subprocess.run(process_cmd, shell=True, check=True)
        
        return output_path

    # process all bed chunks in parallel
    with ThreadPoolExecutor() as executor:
        future_to_file = {executor.submit(process_bed_chunk, file): file for file in temp_bed_files}
        processed_files = []
        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                result = future.result()
                processed_files.append(result)
            except Exception as exc:
                logging.error(f'{file} generated an exception: {exc}')

    # merge processed files
    with tempfile.NamedTemporaryFile(mode='w+t', dir=output_dir, suffix='_merged_3prime_15nt.bed', delete=False) as merged_file:
        merged_path = merged_file.name
    merge_cmd = f"cat {' '.join(processed_files)} | sort -k1,1 -k2,2n > {merged_path}"
    subprocess.run(merge_cmd, shell=True, check=True)

    # intersect with feature files
    results = {}
    for feature, feature_file in [("exon", exon_bed_file), ("intron", intron_bed_file), ("dog", dog_bed_file)]:
        with tempfile.NamedTemporaryFile(mode='w+t', dir=output_dir, suffix=f'_intersect_{feature}.bed', delete=False) as intersect_file:
            intersect_path = intersect_file.name
        intersect_cmd = f"bedtools intersect -a {merged_path} -b {feature_file} -wo -s -sorted -g {genome_file} > {intersect_path}"
        subprocess.run(intersect_cmd, shell=True, check=True)
        
        with open(intersect_path, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                read_id = fields[3]
                feature_info = fields[-3]  # e.g. "ENSG00000136997_region_DOG" or "ENSG00000136997_exon_1"
                # print(f"Feature info: {feature_info}")
                parts = feature_info.split('_')
                if len(parts) >= 3:
                    gene_id = '_'.join(parts[:-2])  
                    feature = '_'.join(parts[-2:])  # feature, e.g., "region_DOG"
                else:
                    gene_id = feature_info  # use whole string if splitting is not applicable
                    feature = "Unknown"  # default feature if splitting does not work

                if read_id not in gene_id_map or gene_id_map[read_id] != gene_id:
                    continue

                overlap = int(fields[-1])
                if read_id not in results:
                    results[read_id] = {}
                if feature not in results[read_id]:
                    results[read_id][feature] = 0
                results[read_id][feature] += overlap


    # retain the best feature per read_id based on maximum overlap
    final_data = []
    for read_id, features in results.items():
        best_feature = None
        max_overlap = -1
        for feature, overlap in features.items():
            if overlap > max_overlap:
                best_feature = (gene_id_map[read_id], feature)
                max_overlap = overlap
        if best_feature:
            final_data.append({'read_id': read_id, 'gene_id': best_feature[0], 'feature': best_feature[1]})

    # create final df 
    df = pd.DataFrame(final_data)
    if not df.empty:
        df.set_index('read_id', inplace=True)  # Set 'read_id' as index if df is not empty
        logging.info("DataFrame created and indexed by read_id.")
    else:
        logging.warning("No valid data processed into DataFrame; DataFrame is empty.")
        
    # print(df.head())

    return df