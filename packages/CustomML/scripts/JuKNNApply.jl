#name: JuKNNApply
#meta.mlname: JuKNN
#meta.mlrole: apply
#description: Custom Julia Apply func for KNN using ScikitLearn
#language: julia
#input: blob model
#input: dataframe df
#input: string namesKeys [original features' names]
#input: string namesValues [new features' names]
#output: dataframe data_out

#>>> Load necessary packages
#using Pkg
#Pkg.add("PyCall")
import PyCall
import BSON: @load
using ScikitLearn

#>>> If original(namesKeys) and new(namesValues) passed, map original names to new
namesValues=split(namesValues, ",")
namesKeys=split(namesKeys, ",")
if length(namesKeys) > 0
    featuresNames = names(df) 
    for i in 1:length(namesKeys)
        replace!(featuresNames, namesValues[i]  =>  namesKeys[i]) 
    end
end
testX=Matrix(select(df, Symbol.(featuresNames)))

#>>> Retrieve saved/trained model
@load file=model trained_model

#>>> Predict using trained model
predY = predict(trained_model, testX)
data_out = DataFrame([predY],[:out])