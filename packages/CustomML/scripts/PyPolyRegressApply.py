#name: PyPolyRegressApply
#meta.mlname: PyPolyRegress
#meta.mlrole: apply
#description: Custom Python Apply func for Polynomail Regression
#language: python
#input: blob model
#input: dataframe df
#input: string namesKeys [original features' names]
#input: string namesValues [new features' names]
#output: dataframe data_out

#>>> Load necessary packages
import numpy as np
import pickle
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

#>>> If original(namesKeys) and new(namesValues) passed, map original names to new
# Temporary commented out block, uncomment when regression is set up
#=============================================================
namesKeys = namesKeys.split(",")
namesValues = namesValues.split(",")
if len(namesKeys) > 0:
    featuresNames = list(df)
    for i in range(len(namesKeys)):
        df = df.rename(columns={namesValues[i]: namesKeys[i]})

testX=np.asarray(df)

#>>> Retrieve saved/trained model
trained_model = pickle.load(open(model, 'rb'))

#>>> Prepare Polynomial test data and predict using trained model
degree = 2 #Temporary hardcoded fix. Degree needs to be retrieve from saved model or passed
poly = PolynomialFeatures(degree = degree)
testX_poly = poly.fit_transform(testX)

#>>> Predict using trained model
predY = trained_model.predict(testX_poly)
data_out = pd.DataFrame({'pred': predY})
