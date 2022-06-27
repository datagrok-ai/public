//name: CountSubsequencePython
//language: javascript
//tags: template, demo
//sample: cars.csv
//input: string sequence
//input: string subsequence
//output: int count
count = 0;
for (let i = 0; i < sequence.length; i++) {
    if (sequence.substring(i, i + subsequence.length) == subsequence) {
        count++;
    }
}