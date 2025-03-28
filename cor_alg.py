#!/usr/bin/env python3
"""
Pipeline for computing correlations and selecting non-redundant probes using GPU.

Script steps:
  1. Computes correlations between probes (samples (rows) x probes (columns)) using GPU
  2. Constructs an adjacency list with correlation values.
  3. Remove redundant probes based on the correlation threshold.
  4. Processes individual chromosome csv files and combines them over multiple iterations.
  
The same steps are executed in several iterations:
  - First iteration: Processes individual chromosomes (chr1, chr2, ..., chr22).
  - Second iteration: Combines 22 directories into 11.
  - Third iteration: Combines 11 directories into 6.
  - Fourth iteration: Combines 6 directories into 3.
  - Final iteration: Combines 3 directories into 1.
  
Usage:
  Example: python3 ./cor_alg.py <dir>
  If <dir> is not provided, the default is "./data".
"""

import os
import pandas as pd
import numpy as np
import cupy as cp
import time
import sys

THRESHOLD = 0.7
BLOCK_SIZE = 1000

def compute_correlation_blockwise(data_cpu, threshold=0.7, block_size=1000):
    """
    Compute correlations between probes using GPU.

    This function computes the correlation matrix:
      - Data is centered and transferred to GPU in blocks.
      - Dot products and norms are computed to obtain correlations.
      - (Optional) Only pairs with |correlation| > threshold are stored in an adjacency list.

    Parameters:
      data_cpu (array): Input probes data (samples (row) x probes (columns)).
      threshold (float): Correlation threshold for storing an edge (Optional use).
      block_size (int): Size of the blocks to be processed on the GPU.

    Returns:
      dict: An adjacency dictionary where each key (probe index) maps to a list of tuples (j, correlation_value)
            for probes j.
    """
    n_samples, n_vars = data_cpu.shape

    # Column means and norms (on CPU)
    col_means = data_cpu.mean(axis=0)
    col_norms = np.empty(n_vars, dtype=np.float32)
    for c in range(n_vars):
        centered_col = data_cpu[:, c] - col_means[c]
        col_norms[c] = np.linalg.norm(centered_col).astype(np.float32)

    # Adjacency list for each probe (Initially empt)
    adjacency = {i: [] for i in range(n_vars)}

    num_col_blocks = (n_vars + block_size - 1) // block_size
    total_block_iters = num_col_blocks * (num_col_blocks + 1) // 2
    iteration_count = 0
    start_time = time.time()

    # Iterate over column blocks (i_block and j_block)
    for i_block in range(num_col_blocks):
        i_start = i_block * block_size
        i_end = min(i_start + block_size, n_vars)

        # Block i on GPU
        block_i_cpu = data_cpu[:, i_start:i_end] - col_means[i_start:i_end]
        block_i_gpu = cp.asarray(block_i_cpu, dtype=cp.float32)
        norms_i_gpu = cp.asarray(col_norms[i_start:i_end], dtype=cp.float32)

        for j_block in range(i_block, num_col_blocks):
            j_start = j_block * block_size
            j_end = min(j_start + block_size, n_vars)

            # Block j on GPU
            block_j_cpu = data_cpu[:, j_start:j_end] - col_means[j_start:j_end]
            block_j_gpu = cp.asarray(block_j_cpu, dtype=cp.float32)
            norms_j_gpu = cp.asarray(col_norms[j_start:j_end], dtype=cp.float32)

            # Dot product between block_i and block_j for correlation
            numer = block_i_gpu.T @ block_j_gpu
            denom = cp.outer(norms_i_gpu, norms_j_gpu) + 1e-20

            corr_block = numer / denom
            # Ensure correlation values lie within [-1, 1]
            corr_block = cp.clip(corr_block, -1.0, 1.0)
            corr_block_cpu = corr_block.get()  # transfer results back to CPU

            rows_block = i_end - i_start
            cols_block = j_end - j_start
            for row in range(rows_block):
                i_global = i_start + row
                for col in range(cols_block):
                    j_global = j_start + col
                    if i_global == j_global:
                        continue  # Skip diagonal elements
                    cval = corr_block_cpu[row, col]
                    # if abs(cval) > threshold:
                    adjacency[i_global].append((j_global, cval))
                    adjacency[j_global].append((i_global, cval))

            # Free GPU memory for block j
            del block_j_gpu, norms_j_gpu, block_j_cpu, numer, denom, corr_block
            cp._default_memory_pool.free_all_blocks()

            iteration_count += 1
            elapsed = time.time() - start_time
            avg_time_per_iter = elapsed / iteration_count
            remaining_iters = total_block_iters - iteration_count
            est_time_left = avg_time_per_iter * remaining_iters
            print(f"Processed block pair: (i_block={i_block}, j_block={j_block}) "
                  f"{iteration_count}/{total_block_iters}, "
                  f"Elapsed: {elapsed:.1f}s, ETA: {est_time_left:.1f}s")

        # Free GPU memory for block i
        del block_i_gpu, norms_i_gpu, block_i_cpu
        cp._default_memory_pool.free_all_blocks()

    return adjacency

def remove_redundant_probes(adjacency, probe_names, threshold):
    """
    Remove redundant probes.

    The algorithm performs iterative removal:
      1. For each active probe, compute the sum of absolute correlations with other active probes.
      2. Sort active probes in descending order of this sum.
      3. For each probe in that order, remove all neighbors with |correlation| > threshold.

    Parameters:
      adjacency (dict): Adjacency list with probe indices and their correlated neighbors.
      probe_names (array): Array of probe names corresponding to indices.
      threshold (float): Correlation threshold.

    Returns:
      tuple:
        - final_remaining (list): List of probe names that were retained.
        - df_iteration_log (df): Summary of each iteration.
        - df_deletion_log (df): Information of each probe removal.
    """
    n = len(probe_names)
    
    # Active probes: True is active, False means removed.
    # Initially, all probes are active.
    active = np.ones(n, dtype=bool)

    deletion_log = []
    iteration_log = []

    start_global = time.time()
    iteration_num = 0

    while True:
        iteration_num += 1
        start_iter = time.time()

        # Get indices active probes
        active_indices = np.where(active)[0]
        if len(active_indices) == 0:
            break

        # 1) Sum of absolute correlations for each active probe
        sum_corr = np.zeros(len(active_indices), dtype=np.float32)
        for idx_k, k in enumerate(active_indices):
            sum_corr[idx_k] = sum(abs(cval) for (j, cval) in adjacency[k] if active[j])

        # 2) Sort active probes in descending order by sum of correlations
        order_sub = np.argsort(-sum_corr)
        removed_this_iter = 0

        # 3) Iterate over the sorted probes
        for pos in order_sub:
            i = active_indices[pos]
            # Skip if probe has already been deactivated
            if not active[i]:
                continue  

            base_probe_name = probe_names[i]
            base_sum_corr = sum_corr[pos]

            # 3.1) Remove highly correlated neighbors
            for (j, cval) in adjacency[i]:
                if active[j] and abs(cval) > threshold:
                    removed_this_iter += 1
                    active[j] = False
                    removed_probe_name = probe_names[j]
                    removed_probe_sum_corr = sum(
                        abs(cval_j) 
                        for (jj, cval_j) in adjacency[j]
                        if active[jj]
                    )
                    
                    deletion_log.append({
                        "base_probe": base_probe_name,
                        "base_probe_sum_corr": base_sum_corr,
                        "removed_probe": removed_probe_name,
                        "removed_probe_sum_corr": removed_probe_sum_corr,
                        "correlation": cval
                    })

        duration_iter = time.time() - start_iter
        iteration_log.append({
            "iteration": iteration_num,
            "num_removed": removed_this_iter,
            "remaining": int(active.sum()),
            "duration_s": round(duration_iter, 2)
        })

        # Terminate flag
        if removed_this_iter == 0:
            break

    duration_total = time.time() - start_global

    print(f"Time for remove redundant probes: {duration_total}")

    df_iteration_log = pd.DataFrame(iteration_log)
    df_deletion_log = pd.DataFrame(deletion_log)
    final_remaining = probe_names[active].tolist()

    return final_remaining, df_iteration_log, df_deletion_log

def process_chromosome_directory(csv_input_path, out_dir, threshold=0.7, block_size=1000):
    """
    Process a chromosome csv file and computing correlations, filtering redundant probes,
    and saving the new data.

    The processing steps are:
      1. Read the input csv.
      3. Ensure that the data matrix has samples as rows and probes as columns 
         (transpose if necessary).
      4. Compute the correlation matrix using GPU acceleration.
      5. Apply the greedy algorithm to remove redundant probes.
      6. Save iteration and deletion logs.
      7. Filter the original data to retain only the selected probes.
      8. Transpose the filtered data (to have probes as rows) and save the final CSV.

    Parameters:
      csv_input_path (str): Path to the input CSV file.
      out_dir (str): Directory where output files will be stored.
      threshold (float): Correlation threshold for probe removal.
      block_size (int): Block size for GPU computation.

    Returns:
      list: List of selected (retained) probe names.
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1) Read the csv file
    df = pd.read_csv(csv_input_path, index_col=0)

    # Remove "chr" column from "split_by_chr.R"
    if "chr" in df.columns:
        df = df.drop(columns=["chr"])

    # 2) For correlation probes have to be the columns and samples the rows, so we make the transpose
    df = df.T  
    # print(df.shape)

    # 3) Convert df to numpy array for correlation calculation
    data_numpy = df.values.astype(np.float32)
    probe_names = df.columns.to_numpy()

    # 4) Compute correlations
    adjacency = compute_correlation_blockwise(
        data_cpu=data_numpy,
        threshold=threshold,
        block_size=block_size
    )

    # 5) Remove redundant probes
    selected_probes, df_iteration_log, df_deletion_log = remove_redundant_probes(
        adjacency=adjacency,
        probe_names=probe_names,
        threshold=threshold
    )

    # 6) Saving results
    results_dir = os.path.join(out_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    iteration_summary_path = os.path.join(results_dir, "iteration_summary.csv")
    detailed_log_path = os.path.join(results_dir, "detailed_removals.csv")
    selected_probes_path = os.path.join(results_dir, "selected_probes.csv")

    df_iteration_log.to_csv(iteration_summary_path, index=False)
    df_deletion_log.to_csv(detailed_log_path, index=False)
    pd.DataFrame({"probe": selected_probes}).to_csv(selected_probes_path, index=False)

    # 7) Filter the original df to retain only the selected probes
    df_filtered = df.loc[:, df.columns.intersection(selected_probes)]

    # 8) Transpose again to have probes as rows and samples as columns
    df_out = df_filtered.T

    # 9) Saving final probes csv file
    name_out = os.path.basename(out_dir)
    csv_out_path = os.path.join(out_dir, f"limma_df_{name_out}.csv")
    df_out.to_csv(csv_out_path)

    print(f"[process_chromosome_directory] Processed and saved results in {out_dir}")
    print(f" - Final number of probes: {len(selected_probes)}\n")

    return selected_probes

def combine_chromosomes(chr_dir_1, chr_dir_2, out_csv_path):
    """
    Combine two chromosome directories by merging the selected probes and their data.

    Steps performed:
      1. Read the 'selected_probes.csv' from both directories.
      2. Create the union of selected probes.
      3. Base names
      4. Load the corresponding csv files (format: probes (rows) x samples (columns)).
      5. Keep the union of probes.
      6. Concatenate data by rows.
      7. Remove duplicate probe entries if any.

    Assumptions:
      - Each csv file is in the format (probes (rows) x samples (columns)) with the index as probes.

    Parameters:
      chr_dir_1 (str): Path to the first chromosome directory.
      chr_dir_2 (str): Path to the second chromosome directory.
      out_csv_path (str): Output file path.
    """
    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)

    # 1) Read the  probes lists from each directory.
    sel_1 = pd.read_csv(os.path.join(chr_dir_1, "results", "selected_probes.csv"))
    sel_2 = pd.read_csv(os.path.join(chr_dir_2, "results", "selected_probes.csv"))

    # 2) Union of all probes
    union_probes = pd.unique(pd.concat([sel_1["probe"], sel_2["probe"]]))

    # 3) Base names
    name_1 = os.path.basename(chr_dir_1)
    name_2 = os.path.basename(chr_dir_2)

    # 4) Load the csv files (format: probes (rows) x samples (columns))
    file_1 = os.path.join(chr_dir_1, f"limma_df_{name_1}.csv")
    file_2 = os.path.join(chr_dir_2, f"limma_df_{name_2}.csv")
    df_1 = pd.read_csv(file_1, index_col=0)
    df_2 = pd.read_csv(file_2, index_col=0)

    # Remove 'chr' column from "split_by_chr.R"
    for col in ["chr"]:
        if col in df_1.columns:
            df_1.drop(columns=[col], inplace=True)
        if col in df_2.columns:
            df_2.drop(columns=[col], inplace=True)

    # 5) Filter rows to keep only the union of selected probes
    df_1_filtered = df_1.loc[df_1.index.intersection(union_probes), :]
    df_2_filtered = df_2.loc[df_2.index.intersection(union_probes), :]

    # 6) Concatenate data by rows
    df_combined = pd.concat([df_1_filtered, df_2_filtered], axis=0)

    # Perhaps a probe is on more than one chromosome?
    # 7) Remove duplicate rows in case a probe appears in more than one file
    df_combined = df_combined.loc[~df_combined.index.duplicated(keep="first")]

    df_combined.to_csv(out_csv_path, index=True)
    print(f"[combine_chromosomes] Combined {name_1} + {name_2} into {out_csv_path}")

def combine_three_chromosomes(chr_dir_1, chr_dir_2, chr_dir_3, out_csv_path):
    """
    Combine three chromosome directories similarly to combine_chromosomes.

    This function reads the selected probes from three directories, merges them,
    loads the corresponding csv files (format: probes (rows) x samples (columns)), filters,
    concatenates, removes duplicates, and saves data.

    Parameters:
      chr_dir_1 (str): Path to the first chromosome directory.
      chr_dir_2 (str): Path to the second chromosome directory.
      chr_dir_3 (str): Path to the third chromosome directory.
      out_csv_path (str): Output file path for the combined CSV.
    """
    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)

    # 1) Read the selected probes lists from each directory
    sel_1 = pd.read_csv(os.path.join(chr_dir_1, "results", "selected_probes.csv"))
    sel_2 = pd.read_csv(os.path.join(chr_dir_2, "results", "selected_probes.csv"))
    sel_3 = pd.read_csv(os.path.join(chr_dir_3, "results", "selected_probes.csv"))

    # 2) Union of all probes
    union_probes = pd.unique(pd.concat([sel_1["probe"], sel_2["probe"], sel_3["probe"]]))

    # 3) Base names
    name_1 = os.path.basename(chr_dir_1)
    name_2 = os.path.basename(chr_dir_2)
    name_3 = os.path.basename(chr_dir_3)

    # 4) Load the csv files (format: probes (rows) x samples (columns))
    file_1 = os.path.join(chr_dir_1, f"limma_df_{name_1}.csv")
    file_2 = os.path.join(chr_dir_2, f"limma_df_{name_2}.csv")
    file_3 = os.path.join(chr_dir_3, f"limma_df_{name_3}.csv")
    df_1 = pd.read_csv(file_1, index_col=0)
    df_2 = pd.read_csv(file_2, index_col=0)
    df_3 = pd.read_csv(file_3, index_col=0)

    # Remove 'chr' column from "split_by_chr.R"
    for col in ["chr"]:
        if col in df_1.columns:
            df_1.drop(columns=[col], inplace=True)
        if col in df_2.columns:
            df_2.drop(columns=[col], inplace=True)
        if col in df_3.columns:
            df_3.drop(columns=[col], inplace=True)

    # 5) Filter rows to keep only the union of selected probes
    df_1_filtered = df_1.loc[df_1.index.intersection(union_probes), :]
    df_2_filtered = df_2.loc[df_2.index.intersection(union_probes), :]
    df_3_filtered = df_3.loc[df_3.index.intersection(union_probes), :]

    # 6) Concatenate data by rows
    df_combined = pd.concat([df_1_filtered, df_2_filtered, df_3_filtered], axis=0)

    # Perhaps a probe is on more than one chromosome?
    # 7) Remove duplicate rows in case a probe appears in more than one file
    df_combined = df_combined.loc[~df_combined.index.duplicated(keep="first")]

    df_combined.to_csv(out_csv_path, index=True)
    print(f"[combine_three_chromosomes] Combined {name_1}, {name_2}, {name_3} into {out_csv_path}")

def run_pipeline(base_dir="./data", threshold=THRESHOLD, block_size=BLOCK_SIZE):
    """
    Pipeline steps:
      1. First iteration: Processing individual chromosomes (chr1, chr2, ..., chr22).
      2. Second iteration: Combining 22 directories into 11.
      3. Third iteration: Combining 11 directories into 6.
      4. Fourth iteration: Combining 6 directories into 3.
      5. Final iteration: Combining 3 directories into 1.

    All operations are performed within the specified base directory.

    Parameters:
      base_dir (str): Base directory containing probes data and where results are saved.
      threshold (float): Correlation threshold for probe removal.
      block_size (int): Block size for GPU correlation computation.
    """
    start_time = time.time()

    # --------------------------
    # 1) First iteration: Process individual chromosomes (chr1, chr2, ..., chr22)
    # --------------------------
    first_iter_base = os.path.join(base_dir, "first_iteration")
    iteration_1_dirs = []

    for i in range(1, 23):
        chr_name = f"chr{i}"
        chr_dir = os.path.join(f"{base_dir}/chr", chr_name)
        
        # csv input file for each chromosome, e.g. "./data/chr1/limma_df_chr1.csv"
        input_csv = os.path.join(chr_dir, f"limma_df_{chr_name}.csv")
        
        out_dir = os.path.join(first_iter_base, chr_name)
        print(f"\n=== First iteration: Processing {chr_name} ===")
        
        # Comment the next line to NOT process each chromosome
        process_chromosome_directory(csv_input_path=input_csv, out_dir=out_dir, threshold=threshold, block_size=block_size)
    
    # --------------------------
    # 2) Second iteration: Combine 22 directories into 11
    # --------------------------
    second_iter_base = os.path.join(base_dir, "second_iteration")
    chr_pairs_2 = [
        ("chr1", "chr22"),
        ("chr2", "chr21"),
        ("chr3", "chr20"),
        ("chr4", "chr19"),
        ("chr5", "chr18"),
        ("chr6", "chr17"),
        ("chr7", "chr16"),
        ("chr8", "chr15"),
        ("chr9", "chr14"),
        ("chr10", "chr13"),
        ("chr11", "chr12"),
    ]
    iteration_2_dirs = []
    for (c1, c2) in chr_pairs_2:
        dir_1 = os.path.join(first_iter_base, c1)
        dir_2 = os.path.join(first_iter_base, c2)
        
        out_name = f"{c1}_{c2}"
        out_dir = os.path.join(second_iter_base, out_name)
        os.makedirs(out_dir, exist_ok=True)
        
        combined_csv = os.path.join(out_dir, f"limma_df_{out_name}.csv")
        
        # Comment to NOT combine directories and process the combined csv:
        combine_chromosomes(chr_dir_1=dir_1, chr_dir_2=dir_2, out_csv_path=combined_csv)
        process_chromosome_directory(csv_input_path=combined_csv, out_dir=out_dir, threshold=threshold, block_size=block_size)
        
        iteration_2_dirs.append(out_name)
    
    # --------------------------
    # 3) Third iteration: Combine 11 directories into 6
    # --------------------------
    third_iter_base = os.path.join(base_dir, "third_iteration")
    # Pairs (or single) from iteration 2 directories
    pairs_3 = [
        (iteration_2_dirs[0], iteration_2_dirs[1]),  # e.g., (chr1_chr22, chr2_chr21)
        (iteration_2_dirs[2], iteration_2_dirs[3]),  # e.g., (chr3_chr20, chr4_chr19)
        (iteration_2_dirs[4], iteration_2_dirs[5]),  # e.g., (chr5_chr18, chr6_chr17)
        (iteration_2_dirs[6], iteration_2_dirs[7]),  # e.g., (chr7_chr16, chr8_chr15)
        (iteration_2_dirs[8], iteration_2_dirs[9]),  # e.g., (chr9_chr14, chr10_chr13)
        (iteration_2_dirs[10], None),                # (chr11_chr12) remains unpaired
    ]
    iteration_3_dirs = []
    
    for (d1, d2) in pairs_3:
        if d2 is None:
            # No combination needed; carry forward directly
            iteration_3_dirs.append(d1)
            continue
        
        new_name = f"{d1}_{d2}"
        out_dir = os.path.join(third_iter_base, new_name)
        os.makedirs(out_dir, exist_ok=True)
        
        dir_1 = os.path.join(second_iter_base, d1)
        dir_2 = os.path.join(second_iter_base, d2)
        
        combined_csv = os.path.join(out_dir, f"limma_df_{new_name}.csv")
        
        # Comment to NOT combine and process:
        combine_chromosomes(chr_dir_1=dir_1, chr_dir_2=dir_2, out_csv_path=combined_csv)
        process_chromosome_directory(csv_input_path=combined_csv, out_dir=out_dir, threshold=threshold, block_size=block_size)
        
        iteration_3_dirs.append(new_name)
    
    # --------------------------
    # 4) Fourth iteration: Combine 6 directories into 3
    # --------------------------
    fourth_iter_base = os.path.join(base_dir, "fourth_iteration")
    pairs_4 = [
        (iteration_3_dirs[0], iteration_3_dirs[1]),
        (iteration_3_dirs[2], iteration_3_dirs[3]),
        (iteration_3_dirs[4], iteration_3_dirs[5] if len(iteration_3_dirs) > 5 else None)
    ]
    iteration_4_dirs = []
    
    for index, (d1, d2) in enumerate(pairs_4):
        if d2 is None:
            iteration_4_dirs.append(d1)
            continue
        
        new_name = f"{d1}_{d2}"
        out_dir = os.path.join(fourth_iter_base, new_name)
        os.makedirs(out_dir, exist_ok=True)
        
        dir_1 = os.path.join(third_iter_base, d1)
        dir_2 = os.path.join(third_iter_base, d2) if os.path.exists(os.path.join(third_iter_base, d2)) else os.path.join(second_iter_base, d2)
        print(dir_2)
        
        combined_csv = os.path.join(out_dir, f"limma_df_{new_name}.csv")
        
        # Comment to NOT combine and process:
        combine_chromosomes(chr_dir_1=dir_1, chr_dir_2=dir_2, out_csv_path=combined_csv)
        process_chromosome_directory(csv_input_path=combined_csv, out_dir=out_dir, threshold=threshold, block_size=block_size)
        
        iteration_4_dirs.append(new_name)
    
    # --------------------------
    # 5) Fifth (final) iteration: Combine 3 directories into 1
    # --------------------------
    final_iter_base = os.path.join(base_dir, "final_iteration")
    os.makedirs(final_iter_base, exist_ok=True)

    if len(iteration_4_dirs) == 3:
        d1, d2, d3 = iteration_4_dirs
        new_name_123 = f"{d1}_{d2}_{d3}"
        out_dir_123 = os.path.join(final_iter_base, new_name_123)
        os.makedirs(out_dir_123, exist_ok=True)
        
        combined_123_csv = os.path.join(out_dir_123, f"limma_df_{new_name_123}.csv")
        
        # Comment to NOT combine three directories and process:
        combine_three_chromosomes(
            chr_dir_1=os.path.join(fourth_iter_base, d1),
            chr_dir_2=os.path.join(fourth_iter_base, d2),
            chr_dir_3=os.path.join(fourth_iter_base, d3),
            out_csv_path=combined_123_csv
        )
        process_chromosome_directory(csv_input_path=combined_123_csv, out_dir=out_dir_123, threshold=threshold, block_size=block_size)

        print(f"Pipeline completed! Final directory: {out_dir_123}")
    else:
        print("There are not exactly 3 combinations in the fourth iteration. Please review the logic.")
    
    elapsed = time.time() - start_time
    print(f"Total execution time: {elapsed}")
    print(iteration_2_dirs)

if __name__ == "__main__":
    """
    Base directory is provided as a command argument, it is used;
    otherwise, the default base directory will be "./data".
    """
    default_base_dir = "./data"

    # command argument for base_dir if provided
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        base_dir = default_base_dir

    print(base_dir)

    run_pipeline(base_dir=base_dir, threshold=0.7, block_size=1024)
