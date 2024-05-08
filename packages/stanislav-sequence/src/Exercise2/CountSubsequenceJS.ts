//name: CountSubsequenceJS 
// input: string sequence
// input: string subsequence
// output: int count 
//language: javascript

export function complementWidget(sequence: any, subsequence: any): any {
    let count = (sequence.match(new RegExp(subsequence, 'gi')) || []).length;
    return count;
}