# -----------------------------------------------------------------
# Legacy exploratory analysis by a member of the Wilson Lab.
# Historical record only — not part of the bancpipeline release.
# -----------------------------------------------------------------

import pandas as pd
import numpy as np
import sqlite3
import os
import json
import pickle
import argparse
import re
import cascade_model

def load_sqlite_database(sql_path):
    """
    Load metadata and connectivity data from SQLite database.
    """
    conn = sqlite3.connect(sql_path)
    meta_df = pd.read_sql_query("SELECT * FROM meta", conn)
    try:
        edgelist_df = pd.read_sql_query("SELECT * FROM edgelist_simple", conn)
    except Exception:
        edgelist_df = pd.read_sql_query("SELECT * FROM edgelist", conn)
    conn.close()
    return meta_df, edgelist_df


def filter_edgelist_and_update_meta(edgelist_df, meta_df, synaptic_threshold=25):
    """
    Filter the edgelist and update metadata based on the synaptic threshold.
    """
    grouped_edgelist = edgelist_df.groupby(['pre', 'post'], as_index=False).agg(synaptic_count=('count', 'sum'))
    filtered_pairs = grouped_edgelist[grouped_edgelist['synaptic_count'] > synaptic_threshold][['pre', 'post']]
    filtered_edgelist = edgelist_df.merge(filtered_pairs, on=['pre', 'post'], how='inner')
    nt_mapping = meta_df.set_index('id')['top_nt'].to_dict()
    filtered_edgelist['top_nt'] = filtered_edgelist['pre'].map(nt_mapping)
    remaining_neurons = set(filtered_edgelist['pre']).union(set(filtered_edgelist['post']))
    updated_meta = meta_df[meta_df['id'].isin(remaining_neurons)].copy()
    filtered_out_neurons = set(meta_df['id']) - remaining_neurons
    return filtered_edgelist, updated_meta, filtered_out_neurons

def process_filtered_edgelist_by_nt(filtered_edgelist, invert_nts=None):
    """
    Process the filtered edgelist by adding a new column 'effective_count',
    which adjusts the 'count' values based on specific neurotransmitters in
    the 'pre_top_nt' column.

    Parameters:
    -----------
    filtered_edgelist : DataFrame
        The edgelist DataFrame with columns including 'pre_top_nt' and 'count'.
    invert_nts : list of str, optional
        List of neurotransmitters for which to invert the 'count' sign.
        Default is ['gaba', 'glutamate'].

    Returns:
    --------
    processed_edgelist : DataFrame
        The modified edgelist with a new column 'effective_count'.
    """
    if invert_nts is None:
        invert_nts = ['gaba', 'glutamate']  # Default neurotransmitters to invert

    # Initialize 'effective_count' as a copy of 'count'
    filtered_edgelist['effective_count'] = filtered_edgelist['count']

    # Process each neurotransmitter in the invert list
    for nt in invert_nts:
        filtered_edgelist.loc[
            filtered_edgelist.top_nt.str.contains(nt, case=False, na=False), 
            "effective_count"
        ] *= -1

    return filtered_edgelist

def extract_mask(meta_df, mask_key, mask_keyword):
    """
    Extract a boolean mask based on the key and keyword.
    """
    if mask_key == "id":
        mask = meta_df[mask_key].astype(str).isin(mask_keyword.split(","))
    else:
        mask = meta_df[mask_key].str.contains(mask_keyword, case=False, na=False)
    return mask


def run_cascade_analysis(meta_df, edgelist_df, mask_A, mask_B, max_timesteps, activation_threshold, n_iterations):
    """
    Run the cascade analysis.
    """
    cascade_simulation = cascade_model.SignalCascade(
        activation_threshold=activation_threshold, 
        n_iterations=n_iterations, 
        max_timesteps=max_timesteps
    )
    start_neurons = set(meta_df[mask_A]['id'])
    end_neurons = set(meta_df[mask_B]['id'])
    result_dict = cascade_simulation.run_cascade_v2(
        S_A0=start_neurons,
        S_E=end_neurons,
        edgelist=edgelist_df,
    )
    return result_dict


def save_results_to_pickle(result_dict, file_path):
    """
    Save the result dictionary to a pickle file.
    """
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'wb') as f:
        pickle.dump(result_dict, f)


def main():
    parser = argparse.ArgumentParser(description="Run cascade analysis with batched queries.")
    parser.add_argument("--sql_path", type=str, 
                        default="/n/data1/hms/neurobio/wilson/banc/connectivity/frankenbrain_v.1.1_data.sqlite",
                        help="Path to the SQLite database.")
    parser.add_argument("--json_path", type=str, default='cascade_batch.json', help="Path to the JSON file with batch parameters.")
    parser.add_argument("--output_dir", type=str, default="cascade_results", help="Directory to save results.")
    args = parser.parse_args()

    # Load SQLite database
    meta_df, edgelist_df = load_sqlite_database(args.sql_path)

    # Load batch parameters from JSON
    with open(args.json_path, 'r') as f:
        batch_params = json.load(f)

    for batch in batch_params:
        synaptic_threshold = batch.get("synaptic_threshold", 25)
        mask_A = batch["mask_A"]
        mask_B = batch["mask_B"]
        mask_key_A = batch.get("mask_key_A", "cell_type")
        mask_key_B = batch.get("mask_key_B", "super_class")
        max_timesteps = batch.get("max_timesteps", 4)
        activation_threshold = batch.get("activation_threshold", 0.01)
        n_iterations = batch.get("n_iterations", 100)

        # Filter metadata and edgelist
        filtered_edgelist, filtered_meta, _ = filter_edgelist_and_update_meta(
            edgelist_df, meta_df, synaptic_threshold
        )
        filtered_edgelist = process_filtered_edgelist_by_nt(filtered_edgelist)
        # Extract masks for start and end neurons
        mask_A_filter = extract_mask(filtered_meta, mask_key_A, mask_A)
        mask_B_filter = extract_mask(filtered_meta, mask_key_B, mask_B)

        # Run cascade analysis
        result_dict = run_cascade_analysis(
            filtered_meta, 
            filtered_edgelist, 
            mask_A_filter, 
            mask_B_filter, 
            max_timesteps, 
            activation_threshold, 
            n_iterations
        )

        # Generate output file path
        mask_A_clean = re.sub(r"\W+", "_", mask_A)
        mask_B_clean = re.sub(r"\W+", "_", mask_B)
        output_filename = f"{mask_A_clean}_to_{mask_B_clean}.pkl"
        output_path = os.path.join(args.output_dir, output_filename)

        # Save results
        save_results_to_pickle(result_dict, output_path)
        print(f"Saved results for {mask_A} to {mask_B} at {output_path}")


if __name__ == "__main__":
    main()