import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import { aarGroups } from '../describe';

export class SubstViewer extends DG.JsViewer {
  viewerGrid: DG.Grid | null;
  maxSubstitutions: number;
  activityLimit: number;
  activityColumnName: string;

  constructor() {
    super();

    this.activityColumnName = this.string('activityColumnName');

    this.maxSubstitutions = this.int('maxSubstitutions', 1);
    this.activityLimit = this.float('acitivityLimit', 0.8);

    this.viewerGrid = null;
  }

  onPropertyChanged(property: DG.Property): void {
    this.calcSubstitutions();
  }

  calcSubstitutions() {
    // let maxSubstitutions: number = 1;
    // let activityLimit = 0.2;

    const aarColName = 'AAR';
    let splitedMatrix: string[][];
    let df: DG.DataFrame = this.dataFrame!;
    const col: DG.Column = df.columns.bySemType('alignedSequence');
    // let values: number[] = df.columns.byName('IC50').toList();
    const values = df.getCol(this.activityColumnName).toList().map(x => -Math.log10(x));
    // values = values;
    splitedMatrix = this.split(col);

    let tableValues: { [aar: string]: number[] } = {};
    let tableFull: { [aar: string]: string[] } = {};

    let nRows = splitedMatrix.length;
    let nCols = splitedMatrix[0].length;

    for (let i = 0; i < nRows - 1; i++) {
      for (let j = i + 1; j < nRows; j++) {
        let substCounter = 0;
        let subst1: { [pos: number]: [string, string] } = {};
        let subst2: { [pos: number]: [string, string] } = {};
        let delta = values[i] - values[j];

        for (let k = 0; k < nCols; k++) {
          if (splitedMatrix[i][k] != splitedMatrix[j][k] && Math.abs(delta) >= this.activityLimit) {
            substCounter++;
            subst1[k] = [splitedMatrix[i][k], splitedMatrix[i][k] + " -> " + splitedMatrix[j][k] + "\t" + values[i] + " -> " + values[j]];
            subst2[k] = [splitedMatrix[j][k], splitedMatrix[j][k] + " -> " + splitedMatrix[i][k] + "\t" + values[j] + " -> " + values[i]];
          }
        }

        if (substCounter <= this.maxSubstitutions && substCounter > 0) {

          Object.keys(subst1).forEach((pos) => {
            let aar = subst1[parseInt(pos)][0];
            if (!Object.keys(tableValues).includes(aar)) {
              tableValues[aar] = Array.apply(null, Array(nCols)).map(function () { return 0; });
              tableFull[aar] = Array.apply(null, Array(nCols)).map(function () { return ""; });
              tableFull[aar][parseInt(pos)] += "Substitution\tvalues\n";
            }

            tableValues[aar][parseInt(pos)]++;
            tableFull[aar][parseInt(pos)] += subst1[parseInt(pos)][1] + "\n";
          });
          Object.keys(subst2).forEach((pos) => {
            let aar = subst2[parseInt(pos)][0];
            if (!Object.keys(tableValues).includes(aar)) {
              tableValues[aar] = Array.apply(null, Array(nCols)).map(function () { return 0; });
              tableFull[aar] = Array.apply(null, Array(nCols)).map(function () { return ""; });
              tableFull[aar][parseInt(pos)] += "Substitution\tvalues\n";
            }

            tableValues[aar][parseInt(pos)]++;
            tableFull[aar][parseInt(pos)] += subst2[parseInt(pos)][1] + "\n";
          });
        }
      }
    }

    const cols = [...Array(nCols).keys()].map((v) => DG.Column.int((v + 1).toString()));
    cols.splice(0, 0, DG.Column.string(aarColName));
    let table = DG.DataFrame.fromColumns(cols);
    for (const aar of Object.keys(tableValues)) {
      table.rows.addNew([aar, ...tableValues[aar]]);
    }

    let groupMapping: { [key: string]: string } = {};

    Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);


    this.viewerGrid = table.plot.grid();

    this.viewerGrid.onCellTooltip(
      (gCell, x, y) => {
        if (gCell.cell.value !== DG.INT_NULL) {
          const colName = gCell.tableColumn!.name;
          if (colName !== aarColName) {
            const aar = this.viewerGrid!.table.get(aarColName, gCell.tableRowIndex!);
            const pos: number = parseInt(colName);
            const tooltipText = tableFull[aar][pos - 1];
            ui.tooltip.show(ui.divText(tooltipText ? tooltipText : 'No substitutions'), x, y);
          }
        }
        return true;
      }
    );

    this.render();
  }

  render() {
    $(this.root).empty();
    this.root.appendChild(this.viewerGrid!.root);
  }

  split(peptideColumn: DG.Column, filter: boolean = true): string[][] {
    const splitPeptidesArray: string[][] = [];
    let currentSplitPeptide: string[];
    let modeMonomerCount = 0;
    let currentLength;
    const colLength = peptideColumn.length;

    // splitting data
    const monomerLengths: { [index: string]: number } = {};
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
    let splitColumns: string[][] = Array.from({ length: modeMonomerCount }, (_) => []);
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
    const columnNames = Array.from({ length: modeMonomerCount }, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1}`);
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

    return splitPeptidesArray;
  }
}
