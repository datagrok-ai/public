import * as DG from 'datagrok-api/dg';

/**
 * Split aligned sequence string into separate parts containing amino acid residues.
 *
 * @export
 * @param {DG.Column} peptideColumn Column containing aligned sequences.
 * @param {boolean} [filter=true] Filter out columns with all the same residues.
 * @return {[DG.DataFrame, number[]]} DataFrame containing split sequence and a list of invalid indexes.
 */
export function splitAlignedPeptides(peptideColumn: DG.Column, filter: boolean = true): [DG.DataFrame, number[]] {
  const splitPeptidesArray: string[][] = [];
  let currentSplitPeptide: string[];
  let modeMonomerCount = 0;
  let currentLength;
  const colLength = peptideColumn.length;

  // splitting data
  const monomerLengths: {[index: string]: number} = {};
  for (let i = 0; i < colLength; i++) {
    currentSplitPeptide = peptideColumn.get(i).split('-').map((value: string) => value ? value : '-');
    splitPeptidesArray.push(currentSplitPeptide);
    currentLength = currentSplitPeptide.length;
    monomerLengths[currentLength + ''] =
      monomerLengths[currentLength + ''] ? monomerLengths[currentLength + ''] + 1 : 1;
  }
  //@ts-ignore: what I do here is converting string to number the most effective way I could find. parseInt is slow
  modeMonomerCount = 1 * Object.keys(monomerLengths).reduce((a, b) => monomerLengths[a] > monomerLengths[b] ? a : b);

  // making sure all of the sequences are of the same size
  // and marking invalid sequences
  let nTerminal: string;
  const invalidIndexes: number[] = [];
  let splitColumns: string[][] = Array.from({length: modeMonomerCount}, (_) => []);
  modeMonomerCount--; // minus N-terminal
  for (let i = 0; i < colLength; i++) {
    currentSplitPeptide = splitPeptidesArray[i];
    nTerminal = currentSplitPeptide.pop()!; // it is guaranteed that there will be at least one element
    currentLength = currentSplitPeptide.length;
    if (currentLength !== modeMonomerCount) {
      invalidIndexes.push(i);
    }
    for (let j = 0; j < modeMonomerCount; j++) {
      splitColumns[j].push(j < currentLength ? currentSplitPeptide[j] : '-');
    }
    splitColumns[modeMonomerCount].push(nTerminal);
  }
  modeMonomerCount--; // minus C-terminal

  //create column names list
  const columnNames = Array.from({length: modeMonomerCount}, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1 }`);
  columnNames.splice(0, 0, 'N-terminal');
  columnNames.push('C-terminal');

  // filter out the columns with the same values
  if (filter) {
    splitColumns = splitColumns.filter((positionArray, index) => {
      const isRetained = new Set(positionArray).size > 1;
      if (!isRetained) {
        columnNames.splice(index, 1);
      }
      return isRetained;
    });
  }

  return [
    DG.DataFrame.fromColumns(splitColumns.map((positionArray, index) => {
      return DG.Column.fromList('string', columnNames[index], positionArray);
    })),
    invalidIndexes,
  ];
}
