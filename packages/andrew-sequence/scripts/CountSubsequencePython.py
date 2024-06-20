# description: counts sequence occurrences
# name: CountSubsequencePython
# language: python
# input: string sequence
# input: string subsequence
# output: int count
count = 0
seq_len = len(sequence)
sub_len = len(subsequence)
for i in range(seq_len - sub_len + 1):
	count += subsequence == sequence[i:i + sub_len]
