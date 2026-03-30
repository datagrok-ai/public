# name: CountSubsequencePythonlocal
# language: python
# input: string sequence
# input: string subsequence
# output: int count


count = 0

for i in range(len(sequence) - len(subsequence) + 1):
    if sequence[i:i + len(subsequence)] == subsequence:
        count += 1
