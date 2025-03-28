#name: Sample Sklearn model 
#description: Trains sklearn linear regression model. Compatible with predictive modeling toolkit
#language: python
#input: dataframe data
#input: string predict_column
#output: blob model [The linear regression model, pickle-serialized]
#meta.queueName: python_docker

from sklearn.linear_model import LinearRegression
import pickle

X = data.drop(predict_column, axis=1)
y = data[predict_column]

lr = LinearRegression()
lr.fit(X, y)
model = pickle.dumps(lr)