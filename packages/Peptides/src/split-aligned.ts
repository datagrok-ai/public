import * as DG from 'datagrok-api/dg';

export function splitAlignedPeptides(peptideColumn: DG.Column) {
  const splitPeptidesArray: string[][] = [];
  let isFirstRun = true;
  let splitted;

  for (const peptideStr of peptideColumn.toList()) {
    splitted = peptideStr.split('-').slice(1, -1);

    if (isFirstRun) {
      for (let i = 0; i < splitted.length; i++) {
        splitPeptidesArray.push([]);
      }
      isFirstRun = false;
    }

    splitted.forEach((value: string, index: number) => {
      splitPeptidesArray[index].push(value === '' ? '-' : value);
    });
  }

  const columnsArray = splitPeptidesArray.map((v, i) => {
    const col = DG.Column.fromList('string', `a${i+1 < 10 ? 0 : ''}${i+1}`, v); //TODO: Remove 'a'
    col.semType = 'aminoAcids';
    return col;
  });

  return DG.DataFrame.fromColumns(columnsArray);
}
