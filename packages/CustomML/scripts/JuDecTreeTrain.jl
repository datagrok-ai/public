#name: JuDecTreeTrain
#meta.mlname: JuDecTree
#meta.mlrole: train
#description: Custom Julia Train func for Decision Tree
#language: julia
#input: dataframe df
#input: string predict_column
#input: int max_depth=2 {category: Parameters}
#output: blob model

#>>> Import necessary packages
using Pkg
Pkg.add("DecisionTree")
using DecisionTree
Pkg.add("BSON")
using BSON: @save

#>>> Extract/Prepare train features and target variable
trainX=Matrix(select(df, Not(Symbol(predict_column))))
trainY=Matrix(select(df, Symbol(predict_column)))
trainY=vec(trainY)

#>>> Build/train model
DTCmodel = DecisionTreeClassifier(max_depth=max_depth)
trained_model = fit!(DTCmodel, trainX, trainY)

#>>> Save trained model
@save file=model trained_model