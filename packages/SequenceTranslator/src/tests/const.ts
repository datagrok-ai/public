type Dict = {[key: string]: string};

export const formatsToHelm: {[key: string]: Dict} = {
  'Axolabs': {
    'UfAfsCfsGfuacg': 'RNA1{[fR](U)p.[fR](A)[sp].[fR](C)[sp].[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)}$$$$'
  },
  'BioSpring': {
    'AT*GC*123456789': 'RNA1{r(A)p.r(T)[sp].r(G)p.r(C)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)p.d([m5C])}$$$$'
  },
  'Mermade12': {
    'hefglijkLIJKHEFG': 'RNA1{[25r](U)[sp].[25r](A)[sp].[25r](C)[sp].[25r](G)[sp].[fR](U)[sp].[fR](A)[sp].[fR](C)[sp].[fR](G)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)}$$$$'
  }
};

export const helmToNucleotides: Dict = {
  'RNA1{[fR](U)p.[fR](A)[sp].[fR](C)[sp].[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)}$$$$': 'UACGUACG',

  'RNA1{r(A)p.r(T)[sp].r(G)p.r(C)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)p.d([m5C])}$$$$': 'ATGCUACGUACGC',

  'RNA1{[25r](U)[sp].[25r](A)[sp].[25r](C)[sp].[25r](G)[sp].[fR](U)[sp].[fR](A)[sp].[fR](C)[sp].[fR](G)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)}$$$$': 'UACGUACGUACGUACG'
};

export const helmToMolfile: Dict = {
};
