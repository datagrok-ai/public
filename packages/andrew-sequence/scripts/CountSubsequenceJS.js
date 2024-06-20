//description: counts sequence occurrences
//name: CountSubsequenceJS
//language: javascript
//input: string sequence
//input: string subsequence
//output: int count
let count = 0;
const seqLen = sequence.length;
const subLen = subsequence.length;

for (let i = 0; i <= seqLen - subLen; ++i)
  count += sequence.slice(i, i + subLen) === subsequence;
