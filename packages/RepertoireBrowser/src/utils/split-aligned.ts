import * as DG from 'datagrok-api/dg';

export function splitAlignedPeptides(peptideColumn: DG.Column) {
  let splitPeptidesArray: string[][] = [];
  let isFirstRun = true;
  let splitted: string[];

  for (const peptideStr of peptideColumn.toList()) {
    splitted = peptideStr.split('-');

    if (isFirstRun) {
      for (let i = 0; i < splitted.length; i++)
        splitPeptidesArray.push([]);

      isFirstRun = false;
    }

    splitted.forEach((value, index) => {
      splitPeptidesArray[index].push(value === '' ? '-' : value);
    });
  }

  //create column names list
  let columnNames = ['N'];
  columnNames = columnNames.concat(splitPeptidesArray.map((_, index) => `${index + 1 < 10 ? 0 : ''}${index +1 }`));
  columnNames.push('C');

  // filter out the columns with the same values
  splitPeptidesArray = splitPeptidesArray.filter((positionArray, index) => {
    const isRetained = new Set(positionArray).size > 1;
    if (!isRetained)
      columnNames.splice(index, 1);

    return isRetained;
  });

  const columnsArray = splitPeptidesArray.map((positionArray, index) => {
    return DG.Column.fromList('string', columnNames[index], positionArray);
  });

  return DG.DataFrame.fromColumns(columnsArray);
}
