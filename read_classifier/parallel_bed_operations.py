from tempfile import NamedTemporaryFile

from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
import concurrent.futures
import subprocess
import csv 

# parent function to split bed file 
def split_bed_file(bed_file, output_dir, num_files):
    # Read file
    with open(bed_file) as f:
        reader = csv.reader(f, delimiter='\t')
        lines = list(reader)

    # Determine chunk size
    chunk_size = len(lines) // num_files
    temp_bed_files = []

    # Create a ThreadPool and submit jobs
    with concurrent.futures.ProcessPoolExecutor(max_workers=48) as executor:
        futures = []
        for i in range(num_files):
            start = i * chunk_size
            end = (i + 1) * chunk_size if i != num_files - 1 else None
            futures.append(executor.submit(write_sorted_chunk, lines[start:end], output_dir, i))

        # Collect the results as they become available
        for future in concurrent.futures.as_completed(futures):
            temp_bed_files.append(future.result())

    return temp_bed_files

# split bedfile into chunks
# daughter function to write sorted chunks 
def write_sorted_chunk(chunk, output_dir, index):
    
    # define the sorted file path
    sorted_file = NamedTemporaryFile(dir=output_dir, delete=False, suffix=f"_sorted_{index}.bed").name

    with open(sorted_file, 'w', newline='') as out:
    
        sort_process = subprocess.Popen(
            ['sort', '-k1,1', '-k2,2n', '--parallel=48', '--buffer-size=80G'],
            stdin=subprocess.PIPE, stdout=out, universal_newlines=True
        )

        for row in chunk:
            print(*row, sep='\t', file=sort_process.stdin)

        # close the sort_process stdin 
        sort_process.stdin.close()

        # wait for the sort process to finish
        sort_process.wait()

    return sorted_file
