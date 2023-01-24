---
title: "Aminoacids WebLogo"
---

Web Logo is used to visualize a graphical representation of multiple sequence alignment (amino acids or
nucleotides). Each logo consists of stacks of symbols, one stack for each position in the sequence.
The overall height of the stack indicates the sequence conservation at that position,
while the height of symbols within the stack indicates the relative frequency of each amino at that position.

In general, a sequence logo provides a richer and more precise description of, for example, a binding site,
than would a consensus sequence.

You must specify the tag `semType` with value `AminoacidsMultipleAlignment` or
`NucleotidesMultipleAlignment` for the data column with multiple alignment sequences, it is mandatory to
select the palette for monomers' colors.

You can customize the look of the viewer with properties. Properties `startPosition` and `endPosition`)
allow to display multiple alignment partially. If property  `startPosition` (`endPosition`)
is not specified, then the Logo will be plotted from the first (till the last) position of sequences.

## General

|             |              |
|-------------|--------------|
| Right click | Context menu |

## Properties

| Property name        | Default | Description                                          |
|----------------------|---------|------------------------------------------------------|
| positionWidth        | 16      | Width of one position stack [px]                     |
| minHeight            | 50      | Minimum height of Logo [px]                          |
| maxHeight            | 100     | Maximum height of Logo [px]                          |
| considerNullSequence | false   | Should logo consider null sequences of data          |
| sequenceColumnName   | null    | source of multiple alignment sequences (column name) |
| startPositionName    | null    | name of the first position to display Logo partially |
| endPositionName      | null    | name of the last position to display Logo partially  |

![Web Logo](web-logo-properties.gif "Web Logo")

## See also

* [Table view](../../datagrok/table-view.md)
* [Viewers](viewers.md)
