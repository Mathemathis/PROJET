import pandas as pd
import matplotlib.pyplot as plt


def influenceEta():
    etas = ["0.01", "0.001", "0.0001", "1.0e-5", "1.0e-6", "1.0e-7"]
    file1_path = "results/heurDich/heurVois" + etas[0] + ".csv" 
    merged_df = pd.read_csv(file1_path)

    for eta in etas[1:]:
        file2_path = "results/heurDich/heurVois" + eta + ".csv" 

        df2 = pd.read_csv(file2_path)

        merged_df = pd.merge(merged_df, df2, on='instance', how='inner', suffixes=("", eta))
    
    print(merged_df.columns)


    return(merged_df)


def boxplot(df):
    plt.figure(figsize=(10, 6))  # Adjust the figure size if needed

    boxplot = df.boxplot(column=['nb_iter', 'nb_iter0.001', 'nb_iter0.0001', 'nb_iter1.0e-5', 'nb_iter1.0e-6', 'nb_iter1.0e-7'],
        showfliers=True, grid=False, patch_artist=True, boxprops=dict(facecolor='lightblue', color='black'))
    boxplot.set_xticklabels(['1e-2', '1e-3', '1e-4', '1e-5', '1e-6', '1e-7'])
    plt.legend()
    plt.title("Variation du nombre d'itérations en fonction de $\eta$")
    plt.xlabel('Valeur de $\eta$')
    plt.ylabel("Nombre d'itérations")
    plt.show()

def countiter(df):
    for eta in ["", "0.001", "0.0001", "1.0e-5", "1.0e-6", "1.0e-7"]:
        print(f"Nombre d'itérations pour eta={eta}: {round(df['nb_iter' + eta].mean(),3)}")

def time(df):
    for eta in ["", "0.001", "0.0001", "1.0e-5", "1.0e-6", "1.0e-7"]:
        print(f"resultat moyen pour eta={eta}: {round(df['ubound' + eta].mean(), 3)}")


def obj(df):
    previous = ""
    for eta in ["0.001", "0.0001", "1.0e-5", "1.0e-6", "1.0e-7"]:
        print("eta = ", eta)
        for index, row in merged_df.iterrows():
            if row['ubound'+previous] != row['ubound'+eta]:
                print(f"Values in ubounds columns are different for instance '{row['instance']}': {row['ubound'+previous]} and {row['ubound'+eta]}")
        previous = eta


def mergevoisins():
    file1_path = "results/heurVois/heurVois3.csv" 
    file2_path = "results/heurVois/heurVois6.csv" 


    df1 = pd.read_csv(file1_path)
    df2 = pd.read_csv(file2_path)


    merged_df = pd.merge(df1, df2, on='instance', how='inner', suffixes=('_3', '_6'))
    return(merged_df)

def analyseVoisinages(df):


    print(f"resultat moyen pour voisinage={3}: {round(df['ubound_3' ].mean(), 3)}")
    res = df['ubound_3' ].mean()
    print(f"resultat moyen pour voisinage={6}: {round(df['ubound_6' ].mean(), 3)}")
    res6 = df['ubound_6' ].mean()
    print(f"Pourcentage d'augmentation: {round(abs(res6-res)/res*100, 2)}%")

    print(f"Temps moyen pour voisinage={3}: {round(df['time_3'].mean(), 3)}")
    res = df['time_3'].mean()
    print(f"Temps moyen pour voisinage={6}: {round(df['time_6'].mean(), 3)}")
    res6 = df['time_6'].mean()
    print(f"Pourcentage d'augmentation: {round(abs(res6-res)/res*100, 2)}%")

    print(f"Nb iterations moyen pour voisinage={3}: {round(df['nb_iter_3'].mean(), 3)}")
    res = df['nb_iter_3'].mean()
    print(f"Nb iterations moyen pour voisinage={6}: {round(df['nb_iter_6'].mean(), 3)}")
    res6 = df['nb_iter_6'].mean()
    print(f"Pourcentage d'augmentation: {round(abs(res6-res)/res*100, 2)}%")

merged_df = mergevoisins()
print(merged_df.columns)
print(merged_df)
analyseVoisinages(merged_df)