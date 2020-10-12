// Computing bitset similarities

let bs1 = DG.BitSet.fromString('00000');
let bs2 = DG.BitSet.fromString('00100');

// should be 1:
console.log(bs1.similarityTo(bs1, "tanimoto"));

// should be 0:
console.log(bs1.similarityTo(bs2, "tanimoto"))