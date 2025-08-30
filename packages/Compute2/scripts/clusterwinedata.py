#name: clusterwinedata
#description: Cluster fetched data
#language: python
#input: dataframe df_wine
#input: int n_pca = 8 { caption: PCA components } [Number of PCA components to cluster]
#input: int n_clusters = 4 {caption: Clusters} [Number of clusters for K-means]
#output: dataframe df_clustering { viewer: scatterPlot(x:"PCA1", y:"PCA2", color:"cluster")}

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

# --- Step 1: Preprocessing ---
features = df_wine.drop(labels=["quality", "wine_type"], axis=1)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(features)

# --- Step 2: PCA for visualization ---
pca = PCA(n_components=n_pca)
X_pca = pca.fit_transform(X_scaled)

# --- Step 3: K-means clustering ---
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
clusters = kmeans.fit_predict(X_scaled).astype(str)

# --- Step 4: Plotting ---
df_clustering = pd.DataFrame({
    "PCA1": X_pca[:, 0],
    "PCA2": X_pca[:, 1],
    "cluster": clusters,
    "wine_type": df_wine["wine_type"]
}, index=df_wine.index)
