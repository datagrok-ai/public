type Dict = {[key: string]: string};

export const formatsToHelm: {[key: string]: Dict} = {
  'Axolabs': {
    'UfAfsCfsGfuacg': 'RNA1{[fR](U)p.[fR](A)[sp].[fR](C)[sp].[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)}$$$$'
  },
  'BioSpring': {
    'AT*GC*123456789': 'RNA1{r(A)p.r(T)[sp].r(G)p.r(C)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)p.d([m5C])}$$$$'
  },
  'Mermade12': {
    'hefglijkLIJKHEFG': 'RNA1{[25r](U)[s.[25r](A)[s.[25r](C)[s.[25r](G)[s.[fR](U)[s.[fR](A)[s.[fR](C)[s.[fR](G)[s.[fR](U)[fR](A)[fR](C)[fR](G)[25r](U)[25r](A)[25r](C)[25r](G)}$$$$'
  }
}

export const helmToNucleotides: Dict = {
  'RNA1{r(A)p.r(T)[sp].r(G)p.r(C)[sp].[fR](U)p.[fR](A)p.[fR](C)p.[fR](G)p.[25r](U)p.[25r](A)p.[25r](C)p.[25r](G)p.d([m5C])}$$$$': 'ATGCUACGUACGC'
}

export const helmToMolfile: Dict = {
}
