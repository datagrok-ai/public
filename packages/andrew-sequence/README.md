# andrew-sequence

`andrew-sequence` is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform (here's a [link](https://dev.datagrok.ai/browse/andrewsequence)).

## Developer exercises

### Exercise 1: Semantic types

+ `complement` & `complementWidget` ([package.ts](./src/package.ts)): the first one takes nucleotides string representation and change each character to the complementary one: A <=> T, G <=> C. The seconds function returns a special info tab widget containing a selected nucleotide data
+ `complement tests` ([sequence-tests.ts](./src/tests/sequence-tests.ts)): cover the `complement` fucntion with a unit test checking dataframes appending and special calculation.
+ `detectNucleotides` detector ([detectors.js](./detectors.js)): defines whether a passed column consists of nucleotide values (strings validation). Returns `dna_nucleotide` semantic type if passed

### Exercise 2-3: Scripting and functions & Composing functions

+ [scripts](./scripts/) dir contains invokable Datagrok scripts in Python & JavaScript:
  + [CountSubsequenceJS](./scripts/CountSubsequenceJS.js) (client function) and [CountSubsequencePython](./scripts/CountSubsequencePython.py) (server function) count sequence occurrences for a passed sequence and by a passed sunsequence
  + [CountSubsequenceJSDataframe](./scripts/CountSubsequenceJSDataframe.js) (client function) and [CountSubsequencePythonDataframe](./scripts/CountSubsequencePythonDataframe.py) (server function) count sequence occurrences for eanch frame values and by a passed sunsequence returning a new column DataFrame with a calculated result column (modifying dataframes)
+ There're also functions in the [main](./src/package.ts) package file to call the scripts using JS API:
  + `callCountSubsequencePythonScript` calls `AndrewSequence:CountSubsequencePython` scripts wrapping its logic (includes some lightweight unit tests)
  + `callCountSubsequenceTableAugmentScript` calls `AndrewSequence:CountSubsequencePythonDataframe` scripts to modify a dataframe (works the same as `CountSubsequencePythonDataframe`)

### Exercise 4: Querying databases

+ [angolovko](./connections/agolovko.json) connection config contains data to handle db access
+ created [queryscript](./queries/queries.sql) `ordersByCountry` can be used to select data from **public.orders** in the **Northwind** storage
+ [package](./src/package.ts) function `getOrders` executes the query script taking `country` valuie as an agrument

### Exercise 5: Reading files

+ [files](./files/) dir contains stored datasets (and other extensions)
+ `openTableFromDemo` function (from [utis](./src/file-units.ts)) opens a dataframe by a passed argument `filepath` using `grok.data.getDemoTable` and adds it to a table view
+ `openTablefromDataFiles` function (from [utis](./src/file-units.ts)) opens a dataframe by a passed argument `filepath` using `grok.data.files.openTable` and adds it to a table view
+ `openTableWithEval` function (from [utis](./src/file-units.ts)) opens a dataframe by a passed argument `filepath` using `grok.functions.eval('OpenServerFile(...)')` and adds it to a table view
+ `addTables` pacakge function works with files distributed with the package, adds all tables (`.csv` files) from the `files` folder to the workspace

### Exercise 6: Creating a scripting viewer

1. A scripting viewer to display a `WEIGHT(HEIGHT)` scatter plot with a color corresponding to `AGE` (from the `Data | Files | Demo Files | demog.csv` dataset):

```python
#name: WEIGHT(HEIGHT) Relation
#language: python
#tags: demo, viewers, hide-suggestions
#input: dataframe t
#input: column xColumnName {type:numerical}
#input: column yColumnName {type:numerical}
#input: column colorColumnName {type:numerical}
#input: column sexColumnName {type:string}
#output: graphics

import numpy as np
import matplotlib.pyplot as plt

color = t[colorColumnName].values
cmap = plt.cm.Spectral
norm = plt.Normalize(vmin=min(color), vmax=max(color))

plt.scatter(t[[xColumnName]], y=t[[yColumnName]], color=cmap(norm(color)), alpha=0.5)
plt.xlabel(xColumnName)
plt.ylabel(yColumnName)
plt.legend()
plt.show()
```

2. A scripting viewer to display a `Frequency of Amino Acid Triplets in Sequences` histogram with all amino acid triplets occurred within all of the sequences in a passed `Sequence` column parameter (from the `Data | Files | Demo Files | bio | sars-cov-2.csv` dataset):

```python
#name: Frequency of Amino Acid Triplets in Sequences
#language: python
#input: dataframe t
#input: column sequenceColumnName {semType: dna_nucleotide}
#output: graphics

from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt

t[sequenceColumnName] = t[sequenceColumnName].str.replace(r'\s+', '', regex=True)

def generate_triplets(sequence):
  return [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]

def generate_triplets_all_frames(sequence):
  triplets = []
  for frame in range(3):
    triplets.extend([sequence[i:i+3] for i in range(frame, len(sequence) - 2, 3)])
  return triplets

triplet_counts = Counter()
for seq in t[sequenceColumnName]:
  triplets = generate_triplets(seq)
  triplet_counts.update(triplets)

plt.figure(figsize=(14, 8))
plt.bar(triplet_counts.keys(), triplet_counts.values())
plt.xlabel('Amino Acid Triplets')
plt.ylabel('Frequency')
plt.title('Frequency of Amino Acid Triplets in Sequences')
plt.xticks(rotation=90)
plt.show()
```

### Exercise 7: Transforming dataframes

+ `fuzzyJoin` pacakge function joins two passed dataframes and adds subsequences `Counts` column (by passed N argument)
+ `fuzzy join function` ([fuzzy-join-test.ts](./src/tests/fuzzy-join-test.ts)): cover the `fuzzyJoin` fucntion with a unit test checking dataframes appending and subsequences calculation.

### Exercise 8: Custom cell renderers

+ `NucleotideBoxCellRenderer` package class renders each sequence nucleotide displaying chars with defined cell padding, offset and a color specified for each nucleotide (`render` method)
+ `nucleotideBoxCellRenderer` package factory function of the `NucleotideBoxCellRenderer` incance for each `dna_nucleotide` row value

### Exercise 9: Creating an info panel with a REST web service

+ `detectNucleotides` detector ([detectors.js](./detectors.js)): defines whether a passed column consists of ENA identifier. Returns `EnaID` semantic type if passed
+ `enaSequence` package function fetches a sequence for the potential corresponding ENA ID in fasta format using **ENA REST API**. Returns a widget containing the request id and a  sequence from a response

### Exercise 10: Enhancing Datagrok with dialog-based functions

+ `_fetchENASequence` package function takes as parameters a `query`, `limit` and `offset` returing a dataframe with two string columns `ID` and `Sequence`. The result is a data from the EMBL database
+ `formENADataTable` package function constructs a dialog form where a `query`, `limit` and `offset` can be configured by a user input to show a fetched table (calling `_fetchENASequence` ) either in a preview mode or table view
