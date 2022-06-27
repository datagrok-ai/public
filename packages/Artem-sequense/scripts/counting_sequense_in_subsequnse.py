#name: CountSubsequencePython
#language: python
#tags: template, demo
#sample: cars.csv
#input: string sequence
#input: string subsequence
#output: int count
count = 0
for i in range(len(sequence)):
    print(sequence[i:i+len(subsequence)])
    if sequence[i:(i+len(subsequence))] == subsequence:
        count += 1
