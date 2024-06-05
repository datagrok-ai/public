# andrew-sequence

`andrew-sequence` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform (here's a [link](https://dev.datagrok.ai/browse/andrewsequence)).

## Developer exercises

### Exercise 1: Semantic types

+ `complement` & `complementWidget` ([package.ts](./src/package.ts)): the first one takes nucleotides string representation and change each character to the complementary one: A <=> T, G <=> C. The seconds function returns a special info tab widget containing a selected nucleotide data.
+ `complement test` ([sequence-tests.ts](./src/tests/sequence-tests.ts)): covers the `complement` fucntion with a unit test checking complemention validity.
+ `detectNucleotides` ([detectors.js](./detectors.js)): efines whether a passed column consists of nucleotide values (strings validation). Returns `dna_nucleotide` semantic type if passed.
