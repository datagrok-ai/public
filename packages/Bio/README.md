# Bio

Bio is a bioinformatics support [package](https://datagrok.ai/help/develop/develop#packages) for the
[Datagrok](https://datagrok.ai) platform with an extensive toolset supporting SAR analisys for small molecules
and antibodies.

## Notations

[@datagrok/bio](https://github.com/datagrok-ai/public/tree/master/packages/Bio) can ingest data in multiple file
formats (such as fasta o csv) and multiple notations for natural and modified residues, aligned and non-aligned forms,
nucleotide and amino acid sequences. The sequences are automatically detected and classified, while preserving their
initial notation. Datagrok allows you to convert sequences between different notations as well.

![Notation converter](../../help/uploads/macromolecules/macromolecules-notation-converter-800.gif "Notation converter")

See:

* [detectMacromolecule()](../Bio/detectors.js)
* [class NotationConverter](../../libraries/bio/src/utils/notation-converter.ts)

## Atomic-Level structures from sequences

For linear sequences, the linear form (see the illustration below) of molecules is reproduced. This is useful
for better visual inspection of sequence and duplex comparison. Structure at atomic level could be saved in available
notations.

![Datagrok-generated atom structure for the ATGCATGC sequence](../../help/uploads/macromolecules/macromolecules-7.png "Datagrok-generated atom structure for the ATGCATGC sequence")

You can easily run this feature for any sequence data using the Bio package and accessing it from the top menu.

![Restoring structure atomic level](../../help/uploads/macromolecules/restoreStructures.gif)

See:

* [getMolfilesFromSeq()](./src/utils/atomic-works.ts)

## MSA

For multiple-sequence alignment, Datagrok uses the “kalign” that relies on Wu-Manber string-matching algorithm
[Lassmann, Timo. _Kalign 3: multiple sequence alignment of large data sets._ **Bioinformatics** (2019).pdf](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz795/30314127/btz795.pdf).
“kalign“ is suited for sequences containing only natural monomers. Sequences of a particular column can be analyzed
using MSA algorithm available at the top menu. Aligned sequences can be inspected for base composition
at the position of MSA result.

![MSA and base composition analysis](
../../help/uploads/macromolecules/macromolecules-msa-and-composition-analysis-800.gif "MSA analysis")

See:

* [runKalign()](src/utils/multiple-sequence-alignment.ts)

TODO: MSA with PepSeA

## Splitting to monomers

Splitting to monomers allows splitting aligned sequences in separate monomers.

![Splitting to monomers](../../help/uploads/macromolecules/splitting-to-monomers.gif)

See:

* [splitAlignedSequences()](../../libraries/bio/src/utils/splitter.ts)

## Web Logo

Web Logo visualizes a graphical representation of multiple sequence alignment (amino acids or nucleotides or
modified residues with multi-char labels). Each logo consists of stacks of symbols, one for each position
in the sequence. The overall height of the stack indicates the sequence conservation at that position,
and the symbol height within the stack indicates the relative frequency of each residue at that position.
In general, a sequence logo provides a more detailed and precise description of, for example, a binding site
than would a consensus sequence.
The most helpful feature for exploration analysis with WebLogo in Datagrok is its ability to control selection
on a dataset. Mouse click on a particular residue in a specific position will select rows of the dataset
with sequences containing that residue at that position.

You must specify the tag ```semType``` with value ```Macromolecule``` and tag `alphabet` of choice ('PT', 'DNA', 'RNA')
for the data column with multiple alignment sequences, it is mandatory to select the palette for monomers' colors.

You can customize the look of the viewer with properties. Properties ```startPosition``` and ```endPosition```)
allow to display multiple alignment partially. If property  ```startPosition``` (```endPosition```)
is not specified, then the Logo will be plotted from the first (till the last) position of sequences.

### General

|             |              |
|-------------|--------------|
| Right click | Context menu |

### Properties

| Property name        | Default  | Description                                                                                                             |
|----------------------|----------|-------------------------------------------------------------------------------------------------------------------------|
| positionWidth        | 16       | Width of one position stack [px]                                                                                        |
| minHeight            | 50       | Minimum height of Logo [px]                                                                                             |
| maxHeight            | 100      | Maximum height of Logo [px]                                                                                             |
| considerNullSequence | false    | Should logo consider null seqences of data                                                                              |
| sequenceColumnName   | null     | source of multiple alignment sequences (column name)                                                                    |
| startPositionName    | null     | name of the first position to display Logo partially                                                                    |
| endPositionName      | null     | name of the last position to display Logo partially                                                                     |
| fixWidth             | false    | Plot takes full width required for sequence length                                                                      |
| verticalAlignment    | 'middle' | choices: ['top', 'middle', 'bottom']                                                                                    |
| horizontalAlignment  | 'center' | choices: ['left', 'center', 'right']                                                                                    |
| fitArea              | true     | Should control to be scaled to fit available area for viewer                                                            |
| shrinkEmptyTail      | true     | Shrink sequences' tails empty in filtered sequences                                                                     |
| skipEmptyPositions   | false    | Skip positions containing only gap symbols in all sequences                                                             |
| positionMarginState  | 'auto'   | choices: ['auto', 'enable', 'off'] Margin between positions. auto - enables margins for sequences of multichar monomers |
| positionMargin       | 0 or 4   | 4 - for sequences of multichar monomers, 0 - single char                                                                |
| positionHeight       | '100%'   | choices: ['100%', 'Entropy'] The way to calculate overall monomers stack height at position                             |

![Web Logo](../../help/visualize/viewers/web-logo-properties.gif "Web Logo")

See also:

* [WebLogo](../../libraries/)
* [Viewers](../../help/visualize/viewers.md)
* [Table view](../../help/datagrok/table-view.md)

## Sequence space

Datagrok allows visualizing multidimensional sequence space using a dimensionality reduction approach.
Several distance-based dimensionality reduction algorithms are available, such as UMAP or t-SNE.
The sequences are projected to 2D space closer if they correspond to similar structures, and farther
otherwise. The tool for analyzing molecule collections is called 'Sequence space' and exists in
the Bio package.

To launch the analysis from the top menu, select Bio | Sequence space.

![Sequence space](../../help/uploads/macromolecules/sequence_space.gif)

See:

* [sequenceSpace()](src/utils/sequence-space.ts)

## Sequence activity cliffs

Activity cliffs tool finds pairs of sequences where small changes in the sequence yield significant
changes in activity or any other numerical property. open the tool from a top menu by selecting.
Similarity cutoff and similarity metric are configurable. As in Sequence space, you can select
from different dimensionality reduction algorithms.

To launch the analysis from the top menu, select Bio | Sequence Activity Cliffs.

![Running activity cliffs](../../help/uploads/macromolecules/activity_cliffs_open.gif)

See:

* [getActivityCliffs()](../../libraries/ml/src/viewers/activity-cliffs.ts)
