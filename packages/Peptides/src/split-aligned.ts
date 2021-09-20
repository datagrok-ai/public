import * as DG from 'datagrok-api/dg';

export function splitAlignedPeptides(peptideColumn: DG.Column) {
  const splitPeptidesArray: string[][] = [];
  let flag = true;
  let splitted;

  for (const peptideStr of peptideColumn.toList()) {
    splitted = peptideStr.split('-').slice(1, -1);

    if (flag) {
      for (let i = 0; i < splitted.length; i++) {
        splitPeptidesArray.push([]);
      }
      flag = false;
    }

    splitted.forEach((value: string, index: number) => {
      splitPeptidesArray[index].push(value === '' ? '-' : value);
    });
  }

  const columnsArray = splitPeptidesArray.map((v, i) => {
    return DG.Column.fromList('string', `${i+1 < 10 ? 0 : ''}${i+1}`, v);
  });

  return DG.DataFrame.fromColumns(columnsArray);
}
