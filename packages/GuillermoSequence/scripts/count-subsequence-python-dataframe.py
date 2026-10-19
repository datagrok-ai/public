#name: CountSubsequencePythonDataframe
#language: python
#input: dataframe sequences
#input: column columnName
#input: string subsequence = "acc"
#output: dataframe result
def count(s):
    return sum(1 for i in range(len(s) - len(subsequence) + 1) if s.startswith(subsequence, i))

result = pd.DataFrame({f'N({subsequence})': sequences[columnName].apply(count)})
