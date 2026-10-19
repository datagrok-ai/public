#name: CountSubsequencePython
#language: python
#input: string sequence
#input: string subsequence
#output: int count
#test: CountSubsequencePython("ATGATC", "A") == 2
#test: CountSubsequencePython("fasta: gttctctacc", "acc") == 1
#test: CountSubsequencePython("aaa", "aa") == 2
count = sum(1 for i in range(len(sequence) - len(subsequence) + 1) if sequence.startswith(subsequence, i))
