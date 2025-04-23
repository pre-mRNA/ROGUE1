#!/usr/bin/env python3
import sys
import re
import array
import logging
import argparse
import multiprocessing
import subprocess
import os

# ---------------------------------------
# Tag-Editing Logic (same as your process_read, but for raw SAM lines)
# ---------------------------------------
def process_line(line):
    """
    Edits a single SAM record's MM/ML tags if present, using logic
    equivalent to the original 'process_read()'. Returns the updated
    SAM line as a string. Header lines pass through unmodified.
    """
    line = line.rstrip("\n")
    if not line or line.startswith("@"):
        return line

    fields = line.split('\t')
    if len(fields) < 11:
        return line

    # optional fields start at index 11
    optional_fields = fields[11:]
    mm_idx = -1
    ml_idx = -1

    # locate MM:Z:... and ML:B:C,...
    for i, field in enumerate(optional_fields):
        if field.startswith("MM:Z:"):
            mm_idx = i
        elif field.startswith("ML:B:C,"):
            ml_idx = i

    # if we don't have both, no modification
    if mm_idx < 0 or ml_idx < 0:
        return line

    mm_str = optional_fields[mm_idx][5:]  # after "MM:Z:"
    ml_str = optional_fields[ml_idx][7:]  # after "ML:B:C,"
    mm_parts = mm_str.split(';')
    try:
        ml_values = list(map(int, ml_str.split(',')))
    except ValueError:
        # malformed ML => skip
        return line

    cumulative_index = 0
    a_plus_a_indices = []
    a_plus_a_nums = []

    # replicate your code's logic
    for part in mm_parts:
        if not part:
            continue
        mod_split = part.split(',')
        mod_name = mod_split[0]
        nums = mod_split[1:]
        num_count = len(nums)
        if mod_name.startswith('A+a'):
            a_plus_a_indices.extend(range(cumulative_index, cumulative_index + num_count))
            a_plus_a_nums.extend([int(n) for n in nums])
        cumulative_index += num_count

    # if no relevant indices, no change
    if not a_plus_a_indices:
        return line

    # boundary check
    if cumulative_index > len(ml_values):
        # mismatch in counts => no change
        return line

    a_plus_a_ml_values = [ml_values[i] for i in a_plus_a_indices]
    filtered_a_plus_a_nums = []
    new_ml_values_filtered = []
    cumulative_sum = 0
    removed_count = 0

    # your original filter logic
    for mod, mlval in zip(a_plus_a_nums, a_plus_a_ml_values):
        if mlval == 255:
            if cumulative_sum > 0 or removed_count > 0:
                filtered_a_plus_a_nums.append(mod + cumulative_sum + removed_count)
                cumulative_sum = 0
                removed_count = 0
            else:
                filtered_a_plus_a_nums.append(mod)
            new_ml_values_filtered.append(mlval)
        else:
            cumulative_sum += mod
            removed_count += 1

    # reconstruct
    if filtered_a_plus_a_nums:
        new_mm_val = "A+a?," + ",".join(map(str, filtered_a_plus_a_nums)) + ";"
        new_ml_val_list = new_ml_values_filtered
    else:
        # remove A+a from mm_parts
        new_mm_parts = [p for p in mm_parts if p and not p.startswith('A+a')]
        new_mm_val = ";".join(new_mm_parts)
        # remove indices from ml
        new_ml_val_list = [val for i, val in enumerate(ml_values) if i not in a_plus_a_indices]

    # rebuild optional fields
    # 1) mm
    if new_mm_val:
        optional_fields[mm_idx] = "MM:Z:" + new_mm_val
    else:
        # remove entirely
        optional_fields.pop(mm_idx)
        # if mm_idx < ml_idx, we need to shift ml_idx
        if mm_idx < ml_idx:
            ml_idx -= 1
    # 2) ml
    if new_ml_val_list:
        new_ml_str = "ML:B:C," + ",".join(map(str, new_ml_val_list))
        if ml_idx < len(optional_fields):
            optional_fields[ml_idx] = new_ml_str
        else:
            # just in case it got popped earlier
            optional_fields.append(new_ml_str)
    else:
        if ml_idx < len(optional_fields):
            optional_fields.pop(ml_idx)

    return "\t".join(fields[:11] + optional_fields)

# ---------------------------------------
# Main Script
# ---------------------------------------
def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

def parse_args():
    parser = argparse.ArgumentParser(description="Filter or edit MM/ML tags in a BAM, saving to new BAM.")
    parser.add_argument("-ibam", required=True, help="input BAM file")
    parser.add_argument("-obam", required=True, help="output BAM file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="number of threads to use")
    parser.add_argument("--chunk-size", type=int, default=500000,
                        help="number of lines in each chunk (for parallel processing)")
    parser.add_argument("--log-every", type=int, default=1000000,
                        help="log progress every N lines (default=1,000,000)")
    return parser.parse_args()

def main():
    args = parse_args()
    setup_logging()

    # sanity checks
    if not os.path.exists(args.ibam):
        logging.error(f"Input BAM does not exist: {args.ibam}")
        sys.exit(1)

    logging.info(f"starting parallel MM/ML tag editing with up to {args.threads} threads")
    logging.info(f"input BAM:  {args.ibam}")
    logging.info(f"output BAM: {args.obam}")

    # samtools pipeline commands
    # 1) convert BAM -> SAM on stdout
    samtools_in_cmd = [
        "samtools", "view", "-h",
        "-@", str(args.threads),
        args.ibam
    ]
    # 2) convert SAM -> BAM on stdin, writing to -obam
    samtools_out_cmd = [
        "samtools", "view",
        "-@", str(args.threads),
        "-b",  # produce BAM
        "-o", args.obam
    ]

    # create a pool of workers
    pool = multiprocessing.Pool(args.threads)

    # use a pipeline
    # We'll read text lines from samtools_in_cmd (SAM) in chunks,
    # process them in parallel, and write lines to samtools_out_cmd.
    total_processed = 0

    try:
        with subprocess.Popen(
            samtools_in_cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True
        ) as in_proc, \
        subprocess.Popen(
            samtools_out_cmd, stdin=subprocess.PIPE, bufsize=1, universal_newlines=True
        ) as out_proc:
            
            chunk = []
            for line in in_proc.stdout:
                # pass header lines directly to output (faster than parallel overhead)
                if line.startswith("@"):
                    out_proc.stdin.write(line)
                    continue
                
                chunk.append(line)
                if len(chunk) >= args.chunk_size:
                    # process chunk in parallel
                    results = pool.map(process_line, chunk, chunksize=10000)
                    for res in results:
                        out_proc.stdin.write(res + "\n")
                    total_processed += len(chunk)

                    # progress indicator
                    if total_processed >= args.log_every and (total_processed % args.log_every) < args.chunk_size:
                        logging.info(f"processed {total_processed} SAM records so far...")

                    chunk.clear()

            # leftover
            if chunk:
                results = pool.map(process_line, chunk, chunksize=10000)
                for res in results:
                    out_proc.stdin.write(res + "\n")
                total_processed += len(chunk)

            # close out_proc.stdin so samtools can finish writing
            out_proc.stdin.close()
            # wait for samtools_out to finish
            out_proc.wait()

        # wait for samtools_in
        in_proc.wait()

    except Exception as e:
        logging.error(f"error during pipeline: {str(e)}", exc_info=True)
        pool.close()
        pool.join()
        sys.exit(1)

    pool.close()
    pool.join()

    logging.info(f"finished processing. total SAM records handled = {total_processed}")

    # index the final BAM
    logging.info(f"indexing output BAM {args.obam} with {args.threads} threads")
    index_cmd = [
        "samtools", "index", "-@", str(args.threads), args.obam
    ]
    ret = subprocess.run(index_cmd)
    if ret.returncode != 0:
        logging.error("samtools index failed!")
        sys.exit(ret.returncode)

    logging.info("all done.")

if __name__ == "__main__":
    main()
