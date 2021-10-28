#name: Chemometric
#language: python
#input: double a
#input: double b
#output: dataframe c1
#output: dataframe c2
#output: int d
#output: double e
#tags: model

d = a + b
e = a / b
d1 = {'col1': [1, 2], 'col2': [3, 4]}
d2 = {'col1': [1, 5], 'col2': [3, 6]}
c1 = pd.DataFrame(data = d1)
c2 = pd.DataFrame(data = d2)
