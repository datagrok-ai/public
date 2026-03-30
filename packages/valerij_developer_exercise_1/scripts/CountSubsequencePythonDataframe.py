#name: CountSubsequencePythonDataframe
#description: Counts subsequence occurrences for each sequence in a dataframe
#language: python
#input: dataframe sequences
#input: column columnName
#input: string subsequence = "acc"
#output: dataframe result {action:join(sequences)}

def count_occurrences(sequence, subsequence):
    sequence = str(sequence).lower().replace('fasta:', '').strip()
    subsequence = str(subsequence).lower().strip()

    count = 0
    for i in range(len(sequence) - len(subsequence) + 1):
        if sequence[i:i + len(subsequence)] == subsequence:
            count += 1
    return count

result = pd.DataFrame()
result[f'N({subsequence})'] = sequences[columnName].apply(
    lambda x: count_occurrences(x, subsequence)
)