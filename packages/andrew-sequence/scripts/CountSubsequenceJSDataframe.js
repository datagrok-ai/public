//description: counts sequence occurrences for frame values
//name: CountSubsequenceJSDataframe
//language: javascript
//input: dataframe sequences
//input: column col {semType: dna_nucleotide}
//input: string subsequence = "acc"
//output: dataframe result {action:join(sequences)}
const countSubsequence = (sequence, subsequence) => {
  let count = 0;
  const seqLen = sequence.length;
  const subLen = subsequence.length;
  for (let i = 0; i <= seqLen - subLen; ++i)
    count += subsequence === sequence.slice(i, i + subLen);
  return count;
};

const subsequenceArr = col.toList().map(x => countSubsequence(x, subsequence));
const resultColName = `N(${subsequence})`;
const resultCol = DG.Column.fromList(DG.TYPE.INT, resultColName, subsequenceArr);
const result = DG.DataFrame.fromColumns([resultCol]);
