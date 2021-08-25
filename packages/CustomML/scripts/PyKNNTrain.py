#name: PyKNNTrain
#meta.mlname: PyKNN
#meta.mlrole: train
#description: Custom Python train func for KNN
#language: python
#input: dataframe df
#input: string predict_column
#input: int n_neighbors {category: FirstParm}
#input: string weights=uniform {category: Parameters; choices: ["uniform", "distance"]}
#input: int leaf_size=30 {category: Parameters}
#input: int p=1 {category: Parameters; range:1-2}
#input: string metric=minkowski {category: Parameters; choices: ["euclidean", "manhattan", "chebyshev", "minkowski"]}
#input: string algorithm = auto {category: Parameters; choices: ["auto", "ball_tree", "kd_tree", "brute"]}
#output: blob model

#>>> Import necessary packages
import numpy as np
import pickle
from sklearn.neighbors import KNeighborsClassifier

#>>> Extract/Prepare train features and target variable
trainX = df.loc[ :,df.columns != predict_column]
trainY = np.asarray (df[predict_column])

#>>> Build and train model
trained_model = KNeighborsClassifier(
    n_neighbors = n_neighbors, 
    weights = weights, 
    leaf_size= leaf_size, 
    p = p, 
    metric = metric, 
    algorithm = algorithm
    )
trained_model.fit(trainX, trainY) 

#>>> Save trained model
pickle.dump(trained_model, open(model, 'wb'))