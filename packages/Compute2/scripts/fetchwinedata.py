#name: fetchwinedata
#description: Fetch Wine quality dataset from web
#language: python
#input: string url_red = "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv"
#input: string url_white = "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
#output: dataframe df_wine
#output: int total_lines

# Load datasets
df_red = pd.read_csv(url_red, sep=';')
df_white = pd.read_csv(url_white, sep=';')

# Add a column to distinguish wine types
df_red["wine_type"] = "red"
df_white["wine_type"] = "white"

# Combine both datasets
df_wine = pd.concat([df_red, df_white], ignore_index=True)
df_wine.reset_index(drop=True, inplace=True)
total_lines = df_wine.shape[0]
