//name: CountSubsequenceJS
//description: Counts occurrences of a subsequence in a sequence
//language: javascript
//input: string sequence
//input: string subsequence
//output: int count

count = 0;

for (let i = 0; i <= sequence.length - subsequence.length; i++) {
  if (sequence.substring(i, i + subsequence.length) === subsequence)
    count++;
}