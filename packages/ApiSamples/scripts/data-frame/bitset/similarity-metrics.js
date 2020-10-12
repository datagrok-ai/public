//Finding the value of similarity

let bs1 = DG.BitSet.create(5);

let bs2 = DG.BitSet.create(5);
bs2.set(0, true)

//should be 1:
console.log(bs1.similarityTo(bs1, "tanimoto"));

//should be 0:
console.log(bs1.similarityTo(bs2, "tanimoto"));