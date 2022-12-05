export const CELL_STRUCTURE = {
  DUPLEX: {
    BEFORE_SS: 'SS ',
    BEFORE_AS: '\r\nAS ',
  },
  TRIPLEX_OR_DIMER: {
    BEFORE_SS: 'SS ',
    BEFORE_AS1: '\r\nAS1 ',
    BEFORE_AS2: '\r\nAS2 ',
  },
};

export function parseStrandsFromDuplexCell(s: string): { SS: string, AS: string } {
  const arr = s
    .slice(CELL_STRUCTURE.DUPLEX.BEFORE_SS.length)
    .split(CELL_STRUCTURE.DUPLEX.BEFORE_AS);
  return {SS: arr[0], AS: arr[1]};
}

export function parseStrandsFromTriplexOrDimerCell(s: string): { SS: string, AS1: string, AS2: string } {
  const arr1 = s
    .slice(CELL_STRUCTURE.TRIPLEX_OR_DIMER.BEFORE_SS.length)
    .split(CELL_STRUCTURE.TRIPLEX_OR_DIMER.BEFORE_AS1);
  const arr2 = arr1[1]
    .split(CELL_STRUCTURE.TRIPLEX_OR_DIMER.BEFORE_AS2);
  return {SS: arr1[0], AS1: arr2[0], AS2: arr2[1]};
}
