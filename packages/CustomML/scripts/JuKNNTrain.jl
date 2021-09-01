#name: JuKNNTrain
#meta.mlname: JuKNN
#meta.mlrole: train
#description: Custom Julia Train func for for KNN using ScikitLearn
#language: julia
#input: dataframe df
#input: string predict_column
#input: int n_neighbors {category: Parameters}
#input: string weights=uniform {category: Parameters; choices: ["uniform", "distance"]}
#input: int leaf_size=30 {category: Parameters}
#input: int p=1 {category: Parameters; range:1-2}
#input: string metric=minkowski {category: Parameters; choices: ["euclidean", "manhattan", "chebyshev", "minkowski"]}
#input: string algorithm = auto {category: Parameters; choices: ["auto", "ball_tree", "kd_tree", "brute"]}
#output: blob model

#>>> Import necessary packages
using Pkg
Pkg.add("ScikitLearn")
using ScikitLearn
@sk_import neighbors: KNeighborsClassifier
Pkg.add("BSON")
using BSON: @save

#>>> Extract/Prepare train features and target variable
trainX=Matrix(select(df, Not(Symbol(predict_column))))
trainY=Matrix(select(df, Symbol(predict_column)))

#>>> Build/train model
trained_model = fit!(KNeighborsClassifier(
    n_neighbors = n_neighbors,
    weights = weights, 
    leaf_size= leaf_size, 
    p = p, 
    metric = metric,
    algorithm = algorithm),
    trainX, trainY)

#>>> Save trained model
@save file=model trained_model