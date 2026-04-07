// Computing bitset similarities

let bs1 = DG.BitSet.fromString('00000');
let bs2 = DG.BitSet.fromString('00100');

// both should be 1:
grok.shell.info(bs1.similarityTo(bs1, 'tanimoto'));
grok.shell.info(bs1.similarityTo(bs1));

// should be 0:
grok.shell.info(bs1.similarityTo(bs2, 'tanimoto'));