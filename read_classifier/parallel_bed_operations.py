import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import csv

def split_bed_file(bed_file, output_dir, num_files):

    temp_bed_files = [os.path.join(output_dir, f"temp_chunk_{i}.bed") for i in range(num_files)]
    file_handlers = [open(temp_file, 'w', newline='') for temp_file in temp_bed_files]
    writers = [csv.writer(f, delimiter='\t') for f in file_handlers]

    with open(bed_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for idx, row in enumerate(reader):
            writers[idx % num_files].writerow(row)

    for f in file_handlers:
        f.close()

    sorted_files = []
    with ProcessPoolExecutor(max_workers=104) as executor:  
        futures = {executor.submit(write_sorted_chunk, chunk_file, output_dir, i): i for i, chunk_file in enumerate(temp_bed_files)}
        for future in as_completed(futures):
            sorted_files.append(future.result())

    for temp_file in temp_bed_files:
        os.remove(temp_file)

    return sorted_files

def write_sorted_chunk(chunk_file, output_dir, index):
    sorted_file = os.path.join(output_dir, f"sorted_chunk_{index}.bed")

    with open(sorted_file, 'w') as out:
        sort_process = subprocess.Popen(
            ['sort', '-k1,1', '-k2,2n', '--parallel=4', '--buffer-size=16G'],  
            stdin=open(chunk_file, 'r'),
            stdout=out
        )
        sort_process.wait()

    return sorted_file
