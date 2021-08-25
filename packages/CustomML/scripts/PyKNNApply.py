#name: PyKNNApply
#meta.mlname: PyKNN
#meta.mlrole: apply
#description: Custom Python Apply for KNN
#language: python
#input: blob model
#input: dataframe df
#input: string namesKeys [original features' names]
#input: string namesValues [new features' names]
#output: dataframe data_out

#>>> Load necessary packages
import numpy as np
import pickle

#>>> If original(namesKeys) and new(namesValues) passed, map original names to new
namesKeys = namesKeys.split(",")
namesValues = namesValues.split(",")
if len(namesKeys) > 0:
    featuresNames = list(df)
    for i in range(len(namesKeys)):
        df = df.rename(columns={namesValues[i]: namesKeys[i]})

testX=np.asarray(df)

#>>> Retrieve saved/trained model
trained_model = pickle.load(open(model, 'rb'))

#>>> Predict using trained model
predY = trained_model.predict(testX)
data_out = pd.DataFrame({'pred': predY})
