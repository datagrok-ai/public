#name: PyPolyRegressTrain
#meta.mlname: PyPolyRegress
#meta.mlrole: train
#description: Custom Python train func for Polynomail Regression
#language: python
#input: dataframe df
#input: string predict_column
#input: int degree=2 {category: Parameters}
#output: blob model

#>>> Import necessary packages
import numpy as np
import pickle
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

#>>> Extract/Prepare train features and target variable
trainX = df.loc[ :,df.columns != predict_column]
trainY = np.asarray (df[predict_column])

#>>> Build and train model
poly = PolynomialFeatures(degree = degree)
trainX_poly = poly.fit_transform(trainX)
trained_model = LinearRegression().fit(trainX_poly, trainY) 

#>>> Save trained model
pickle.dump(trained_model, open(model, 'wb'))

