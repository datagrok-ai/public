# description: counts sequence occurrences for frame values
# name: CountSubsequencePythonDataframe
# language: python
# input: dataframe sequences
# input: column columnName {semType: dna_nucleotide}
# input: string subsequence = "acc"
# output: dataframe result {action:join(sequences)}
def count_subsequence(sequence, subsequence):
    count = 0
    seq_len = len(sequence)
    sub_len = len(subsequence)
    for i in range(seq_len - sub_len + 1):
        count += subsequence == sequence[i:i + sub_len]
    return count

result_series = sequences[columnName].apply(lambda x: count_subsequence(x, subsequence))
result = pd.DataFrame({f'N({subsequence})': result_series})
