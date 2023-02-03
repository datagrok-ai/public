// Bitwise AND, OR, XOR, and AND_NOT operations

let bs1 = DG.BitSet.fromString('11100');
let bs2 = DG.BitSet.fromString('00111');

let s = `
  b1: ${bs1.toBinaryString()}
  b2: ${bs2.toBinaryString()}
  and: ${bs1.clone().and(bs2).toBinaryString()}
  or: ${bs1.clone().or(bs2).toBinaryString()}
  xor: ${bs1.clone().xor(bs2).toBinaryString()}
  and not: ${bs1.clone().andNot(bs2).toBinaryString()}
`;

grok.shell.info(s);