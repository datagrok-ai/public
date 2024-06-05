# andrew-sequence

`andrew-sequence` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform (here's a [link](https://dev.datagrok.ai/browse/andrewsequence)).

## Developer exercises

### Exercise 1: Semantic types

+ `complement` & `complementWidget` ([package.ts](./src/package.ts)): the first one takes nucleotides string representation and change each character to the complementary one: A <=> T, G <=> C. The seconds function returns a special info tab widget containing a selected nucleotide data
+ `complement testы` ([sequence-tests.ts](./src/tests/sequence-tests.ts)): cover the `complement` fucntion with unit testы checking complemention validity cases.
+ `detectNucleotides` ([detectors.js](./detectors.js)): efines whether a passed column consists of nucleotide values (strings validation). Returns `dna_nucleotide` semantic type if passed

### Exercise 2-3: Scripting and functions & Composing functions

+ [scripts](./scripts/) dir contains invokable Datagrok scripts in Python & JavaScript:
  + [CountSubsequenceJS](./scripts/CountSubsequenceJS.js) (client function) and [CountSubsequencePython](./scripts/CountSubsequencePython.py) (server function) count sequence occurrences for a passed sequence and by a passed sunsequence
  + [CountSubsequenceJSDataframe](./scripts/CountSubsequenceJSDataframe.js) (client function) and [CountSubsequencePythonDataframe](./scripts/CountSubsequencePythonDataframe.py) (server function) count sequence occurrences for eanch frame values and by a passed sunsequence returning a new column DataFrame with a calculated result column (modifying dataframes)
+ There're also functions in the [main](./src/package.ts) package file to call the scripts using JS API:
  + `callCountSubsequencePythonScript` calls `AndrewSequence:CountSubsequencePython` scripts wrapping its logic (includes some lightweight unit tests)
  + `callCountSubsequenceTableAugmentScript` calls `AndrewSequence:CountSubsequencePythonDataframe` scripts to modify a dataframe (works the same as `CountSubsequencePythonDataframe`)
