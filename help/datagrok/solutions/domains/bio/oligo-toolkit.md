---
title: "Oligo Toolkit"
sidebar_position: 1
---

The OligoToolkit is a collection of tools helping you to work with oligonucleotide sequences 
– primarily for synthesizing DNA/RNA oligonucleotides with various chemical modifications.
An example use-case is producing double-stranded RNAs for the RNA-silencing.

To open the app, on the **Sidebar**, select **Functions** > **Apps** > **Oligo Toolkit**.

![Run OligoToolkit](img/Oligotoolkit-run.gif)

The OligoToolkit 
application contains several independent modules, 
allowing you to do the following:

* Individual and batch modification of the FASTA-encoded DNA/RNA sequences according to the pattern. 
* Conversion of chemically modified sequences between several open and commercial formats: FASTA, HELM, MerMade, BioSpring.
* Producing chemical structures for the natural and modified oligonucleotides.

## Pattern

The Pattern tool applies modification patterns to a set of sequences.
1. You specify the pattern: create it manually or load it from saved patterns. 
   Loading patterns created by other users are also available.
   For the pattern, you can choose the sequence basis (DNA/RNA) and pattern length
   for sense and anti-sense strands.
2. You upload the Excel/CSV file containing unmodified sequences — the format: sequence ID, sense, and antisense strands. 
3. You convert it and receive the modified sequence in HELM or custom synthesizer format.  
4. The conversion utility designed for dimers: By default, you specify patterns for sense and anti-sense strands. 
   To disable the generation of the complementary chain, disable the checkbox below the pattern window.

**Available modifications:**

* Type of the bases: Methylated, with 2’-fluorescent mark, artificial bases.
* Additional groups at the 5’ and 3’ ends of the chain.
* Phospho-tiole connections.
* Nucleotide connection type (2’ instead of 3’).
* Tio-sulfur bond between nucleotides.



## Translator

A simple tool to convert oligonucleotide sequences 
between different formats and generate SDF and SMILES representations. 

You specify one sequence, choose a format, and receive the set of sequences in all other supported formats. 
The Translator application has no bulk load supported. 
This tool is designed to convert only 1-2 sequences. 
The most common use case is checking that the conversion is correct for a particular structure.

![Oligotoolkit-Translator](img/Oligotoolkit-translator.gif)

### Monomer library

Click the book icon next to the caption on the Translator tab to open the Monomer library – 
the collection of all monomers that the Translator uses to generate chemical formulas.
The monomers used by sequencing machines do not follow the monomer definition in the HELM standard. 
So, historically, the OligoToolkit monomer library is a simple dataframe 
containing monomer chemical structures and codes for sequence formats.
To download the library, click the “Down arrow” icon next to the caption and select format. 
Currently, the only way to add new monomers to this library and specify connection points 
is by using the custom Python script from the OligoToolkit package, which is tricky.
This part of the code is being refactored now.

### Structure

A simple tool to convert sequences to SDF. 
It is designed especially for generating siRNA chemical structure 
with the sense strand and one or two anti-sense strands. 
Paste the sequences for sense and anti-sense strands and get the structure in the SDF format.
You can save separate structures for each strand or combined SDF for all strands.
Also, you can choose the direction for strands.

![Oligotoolkit-Structure](img/Oligotoolkit-structure.gif)

## OligoBatch Calculator

Separate specialized applications, calculating the properties of the oligonucleotides 
(for example, obtained from sequencing machine or purchased).
Specify the yield amount: value and units. 
Then, paste oligonucleotide sequences in the FASTA or AxoLabs format.
The tool calculates the table containing oligonucleotide properties: 
length, molecular weight, the amount in different units (mkg, nMols), optical density, and extinction coefficient.
