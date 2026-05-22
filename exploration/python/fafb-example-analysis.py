# -----------------------------------------------------------------
# Legacy exploratory analysis by a member of the Wilson Lab.
# Historical record only — not part of the bancpipeline release.
# -----------------------------------------------------------------

import os
import sqlite3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

# Data sources
fafb_sql = "/Volumes/neurobio/wilsonlab/banc/connectivity/flywire_783_data.sqlite"
hemibrain_sql = "/Volumes/neurobio/wilsonlab/banc/connectivity/hemibrain_v.1.2.1_data.sqlite"
manc_sql = "/Volumes/neurobio/wilsonlab/banc/connectivity/manc_1.2.1_data.sqlite"
banc_sql = "/Volumes/neurobio/wilsonlab/banc/connectivity/banc_data.sqlite"

# Output
fafb_images = "inst/images/fafb/"

# Connect to the database
conn = sqlite3.connect(fafb_sql)

# List all tables in the database
cursor = conn.cursor()
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
print("Tables:", tables)

# Get the meta data, cell types, etc.
fw_meta = pd.read_sql_query("SELECT * FROM meta", conn)
fw_meta = fw_meta.sort_values('top_p', ascending=False).groupby('cell_type').first().reset_index()

# Choose analysis problem
fw_tangential = fw_meta[fw_meta['cell_type'].str.contains('CB.FB') & 
                        fw_meta['cell_class'].isin(['central_complex', 'CX'])]
fw_hdelta = fw_meta[fw_meta['cell_type'].str.contains('hDelta') & 
                    fw_meta['cell_class'].isin(['central_complex', 'CX'])]

tangential_ids = fw_tangential['root_783'].tolist()
hdelta_ids = fw_hdelta['root_783'].tolist()

# Construct the SQL query
query = f"""
SELECT *
FROM edgelist
WHERE count >= 20
  AND pre IN ({','.join(['?']*len(tangential_ids))})
  AND post IN ({','.join(['?']*len(hdelta_ids))})
"""

# Execute the query and collect the results
fw_elist_th = pd.read_sql_query(query, conn, params=tangential_ids + hdelta_ids)

# Merge the edgelist with meta data
fw_elist_merged = fw_elist_th.merge(fw_meta[['root_783', 'cell_type', 'top_nt']], 
                                    left_on='pre', right_on='root_783', 
                                    suffixes=('', '_pre'))
fw_elist_merged = fw_elist_merged.merge(fw_meta[['root_783', 'cell_type', 'top_nt']], 
                                        left_on='post', right_on='root_783', 
                                        suffixes=('_pre', '_post'))
fw_elist_merged = fw_elist_merged.rename(columns={'cell_type_pre': 'pre_cell_type', 
                                                  'cell_type_post': 'post_cell_type', 
                                                  'top_nt_pre': 'pre_top_nt', 
                                                  'top_nt_post': 'post_top_nt'})

# Reshape and collapse the data by cell_type
fw_elist_collapsed = fw_elist_merged.groupby(['pre_cell_type', 'post_cell_type', 
                                              'pre_label', 'post_label']).agg({
    'count': 'mean',
    'pre_top_nt': 'first',
    'post_top_nt': 'first'
}).reset_index()
fw_elist_collapsed['connection_type'] = fw_elist_collapsed['pre_label'] + ' to ' + fw_elist_collapsed['post_label']

def create_heatmap(data, conn_type, filename):
    mat_data = data[data['connection_type'] == conn_type]
    
    # Create matrix
    mat = mat_data.pivot(index='pre_cell_type', columns='post_cell_type', values='count').fillna(0)
    
    # Convert matrix to numeric
    mat = mat.astype(float)
    
    # Format numbers for display
    mat_display = mat.applymap(lambda x: f'{x:.2f}')
    
    # Create annotation data frames
    row_ann = mat_data.groupby('pre_cell_type')['pre_top_nt'].first()
    col_ann = mat_data.groupby('post_cell_type')['post_top_nt'].first()
    
    # Create color palette for top_nt
    nt_colors = {
        'acetylcholine': "#EF7C12",
        'glutamate': "#8FDA04",
        'gaba': "#1BB6AF",
        'serotonin': "#FBBB48",
        'dopamine': "#F4E3C7",
        'octopamine': "#C70E7B",
        'unknown': "grey"
    }
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(mat, annot=mat_display, fmt='', cmap='YlOrRd', ax=ax)
    
    # Add row and column annotations
    ax2 = ax.twinx()
    ax3 = ax.twiny()
    
    for i, nt in enumerate(row_ann):
        ax2.add_patch(plt.Rectangle((0, i), 0.05, 1, color=nt_colors.get(nt, 'grey')))
    
    for i, nt in enumerate(col_ann):
        ax3.add_patch(plt.Rectangle((i, mat.shape[0]), 1, 0.05, color=nt_colors.get(nt, 'grey')))
    
    ax2.set_ylim(ax.get_ylim())
    ax3.set_xlim(ax.get_xlim())
    ax2.axis('off')
    ax3.axis('off')
    
    plt.title(conn_type)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# Create heatmaps for each connection type
for conn_type in fw_elist_collapsed['connection_type'].unique():
    filename = os.path.join(fafb_images, f"heatmap_collapsed_{conn_type.replace(' ', '_')}.png")
    create_heatmap(fw_elist_collapsed, conn_type, filename)

# Close the database connection
conn.close()