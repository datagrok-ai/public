// Getting BitSet from bytes

let arr8 = new Uint8Array(4);
arr8[0] = 1;
arr8[2] = 1;

let bs1 = DG.BitSet.fromBytes(arr8.buffer);

// should be 10000000000000001000000000000000:
console.log(bs1.toBinaryString());