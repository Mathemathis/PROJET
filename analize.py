import pandas as pd

file1_path = "results/heurVois/resultsHeurVois.csv" 
file2_path = "results/heurVois/resultsHeurVois2.csv" 


df1 = pd.read_csv(file1_path)
df2 = pd.read_csv(file2_path)


merged_df = pd.merge(df1, df2, on='instance', how='inner', suffixes=('_1', '_2'))

for index, row in merged_df.iterrows():
    if row['ubound_1'] != row['ubound_2']:
        print(f"Values in ubounds columns are different for instance '{row['instance']}': {row['ubound_1']} and {row['ubound_2']}")