# name: CountSubsequencePythonDataframe
# language: python
# input: dataframe sequences
# input: column columnName
# input: string subsequence = "acc"
# output: dataframe result {action:join(sequences)}

data = pd.DataFrame(sequences)[columnName];
subsequenceCounts = []
for x in data:
  subsequenceCounts.append(x.count(subsequence));

result = pd.DataFrame({'N('+subsequence+')': subsequenceCounts});
print(result)
