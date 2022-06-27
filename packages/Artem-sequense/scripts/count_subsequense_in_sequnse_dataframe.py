#name: CountSubsequencePythonDataframe
#language: python
#input: dataframe sequences
#input: column columnName
#input: string subsequence = "acc"
#output: dataframe result {action:join(sequences)}
sequencesColumn = sequences[columnName]
count = 0
subsequence = "a"
result = pd.DataFrame(columns=['N(s)'])
for i in range(len(sequencesColumn)):
    count = 0
    for j in range(len(sequencesColumn[i])):
        if sequencesColumn[i][j:j+len(subsequence)] == subsequence:
            count += 1
    row_df = pd.DataFrame({'N(s)': [count]})
    result = pd.concat([result, row_df], ignore_index=True)
print(result)

