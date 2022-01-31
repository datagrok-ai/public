import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

// import {aarGroups} from '../describe';
import {setAARRenderer} from '../utils/cell-renderer';
import { SemanticValue } from 'datagrok-api/dg';

export class SubstViewer extends DG.JsViewer {
  viewerGrid: DG.Grid | null;
  maxSubstitutions: number;
  activityLimit: number;
  activityColumnName: string;
  title: string;
  casesGrid: DG.Grid | null;

  constructor() {
    super();

    this.activityColumnName = this.string('activityColumnName');
    this.title = this.string('title', 'Substitution analysis');

    this.maxSubstitutions = this.int('maxSubstitutions', 1);
    this.activityLimit = this.float('activityLimit', 2);

    this.viewerGrid = null;
    this.casesGrid = null;
  }

  onPropertyChanged(property: DG.Property): void {
    this.calcSubstitutions();
  }

  calcSubstitutions() {
    const aarColName = 'AAR';
    const df: DG.DataFrame = this.dataFrame!;
    const col: DG.Column = df.columns.bySemType('alignedSequence');
    // let values: number[] = df.columns.byName('IC50').toList();
    const values = df.getCol(this.activityColumnName).toList();
    // values = values;
    const splitedMatrix = this.split(col);

    const tableValues: { [aar: string]: number[] } = {};
    const tableTooltips: { [aar: string]: {[index: string]: string}[][] } = {};
    const tableCases: { [aar: string]: number[][][] } = {};

    const nRows = splitedMatrix.length;
    const nCols = splitedMatrix[0].length;
    const nColsArray = Array(nCols);

    for (let i = 0; i < nRows - 1; i++) {
      for (let j = i + 1; j < nRows; j++) {
        let substCounter = 0;
        const subst1: { [pos: number]: [string, {[index: string]: string}] } = {};
        const subst2: { [pos: number]: [string, {[index: string]: string}] } = {};
        const delta = values[i] - values[j];

        for (let k = 0; k < nCols; k++) {
          const smik = splitedMatrix[i][k];
          const smjk = splitedMatrix[j][k];
          if (smik != smjk && Math.abs(delta) >= this.activityLimit) {
            const vi = values[i].toFixed(2);
            const vj = values[j].toFixed(2);
            substCounter++;
            subst1[k] = [
              smik,
              {key: `${smik === '-' ? 'Empty' : smik} → ${smjk === '-' ? 'Empty' : smjk}`, value: `${vi} → ${vj}`},
            ];
            subst2[k] = [
              smjk,
              {key: `${smjk === '-' ? 'Empty' : smjk} → ${smik === '-' ? 'Empty' : smik}`, value: `${vj} → ${vi}`},
            ];
          }
        }

        if (substCounter <= this.maxSubstitutions && substCounter > 0) {
          Object.keys(subst1).forEach((pos) => {
            const posInt = parseInt(pos);
            const aar = subst1[posInt][0];
            if (!Object.keys(tableValues).includes(aar)) {
              tableValues[aar] = Array.apply(null, nColsArray).map(function() {
                return DG.INT_NULL;
              });
              tableTooltips[aar] = Array.apply(null, nColsArray).map(function() {
                return [];
              });
              tableCases[aar] = Array.apply(null, nColsArray).map(function() {
                return [];
              });
            }

            tableValues[aar][posInt] = tableValues[aar][posInt] === DG.INT_NULL ? 1 : tableValues[aar][posInt] + 1;
            tableTooltips[aar][posInt] = !tableTooltips[aar][posInt].length ?
              [{key: 'Substitution', value: 'Values'}] : tableTooltips[aar][posInt];
            tableTooltips[aar][posInt].push(subst1[posInt][1]);
            tableCases[aar][posInt].push([i, j, delta]);
          });
          Object.keys(subst2).forEach((pos) => {
            const posInt = parseInt(pos);
            const aar = subst2[posInt][0];
            if (!Object.keys(tableValues).includes(aar)) {
              tableValues[aar] = Array.apply(null, nColsArray).map(function() {
                return DG.INT_NULL;
              });
              tableTooltips[aar] = Array.apply(null, nColsArray).map(function() {
                return [];
              });
              tableCases[aar] = Array.apply(null, nColsArray).map(function() {
                return [];
              });
            }

            tableValues[aar][posInt] = tableValues[aar][posInt] === DG.INT_NULL ? 1 : tableValues[aar][posInt] + 1;
            // tableValues[aar][posInt]++;
            tableTooltips[aar][posInt] = !tableTooltips[aar][posInt].length ?
              [{key: 'Substitution', value: 'Values'}] : tableTooltips[aar][posInt];
            tableTooltips[aar][posInt].push(subst2[posInt][1]);
            tableCases[aar][posInt].push([j, i, -delta]);
          });
        }
      }
    }

    const tableValuesKeys = Object.keys(tableValues);
    const dfLength = tableValuesKeys.length;
    const cols = [...nColsArray.keys()].map((v) => DG.Column.int(v.toString(), dfLength));
    cols.forEach((col: DG.Column) => col.semType = 'Substitution');
    const aarCol = DG.Column.string(aarColName, dfLength);
    cols.splice(0, 1, aarCol);
    const table = DG.DataFrame.fromColumns(cols);

    for (let i = 0; i < dfLength; i++) {
      const aar = tableValuesKeys[i];
      tableValues[aar].splice(0, 1);
      table.rows.setValues(i, [aar, ...tableValues[aar]]);
    }

    // let groupMapping: { [key: string]: string } = {};

    //TODO: enable grouping
    // Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);

    this.viewerGrid = table.plot.grid();

    setAARRenderer(aarCol, this.viewerGrid);

    this.viewerGrid.onCellTooltip(
      (gCell, x, y) => {
        if (gCell.cell.value !== DG.INT_NULL && gCell.tableColumn !== null && gCell.tableRowIndex !== null) {
          const colName = gCell.tableColumn.name;
          if (colName !== aarColName) {
            const aar = this.viewerGrid!.table.get(aarColName, gCell.tableRowIndex);
            const pos = parseInt(colName);
            const tooltipText = tableTooltips[aar][pos].length ?
              DG.HtmlTable.create(
                tableTooltips[aar][pos], (item: {[index: string]: string}, idx: number) => [item.key, item.value],
              ).root : ui.divText('No substitutions');
            ui.tooltip.show(tooltipText, x, y);
          }
        }
        return true;
      },
    );

    for (const col of table.columns.names())
      this.viewerGrid.col(col)!.width = this.viewerGrid.props.rowHeight;


    this.viewerGrid.onCellRender.subscribe((args) => {
      if (args.cell.isRowHeader && args.cell.gridColumn.visible) {
        args.cell.gridColumn.visible = false;
        args.preventDefault();
      }
    });

    this.viewerGrid.props.allowEdit = false;

    table.onCurrentCellChanged.subscribe((_) => {
      if (table.currentCol !== null && table.currentCol.name !== aarColName && table.currentCell.value !== null) {
        const aar = table.get(aarColName, table.currentRowIdx);
        const pos = parseInt(table.currentCol.name);
        const currentCase = tableCases[aar][pos];
        const tempDfLength = currentCase.length;
        const initCol = DG.Column.string('Initial', tempDfLength);
        const subsCol = DG.Column.string('Substituted', tempDfLength);

        const tempDf = DG.DataFrame.fromColumns([
          initCol,
          subsCol,
          DG.Column.float('Difference', tempDfLength),
        ]);

        for (let i = 0; i < tempDfLength; i++) {
          const row = currentCase[i];
          tempDf.rows.setValues(i, [col.get(row[0]), col.get(row[1]), row[2]]);
        }

        initCol.semType = 'alignedSequence';
        initCol.setTag('isAnalysisApplicable', 'false');
        subsCol.semType = 'alignedSequence';
        subsCol.setTag('isAnalysisApplicable', 'false');

        this.casesGrid = tempDf.plot.grid();
        this.casesGrid.props.allowEdit = false;
        grok.shell.o = SemanticValue.fromValueType(tempDf, 'Substitution');
      } else {
        grok.shell.o = SemanticValue.fromValueType(null, 'Substitution');
        this.casesGrid = null;
      }

      this.render();
    });

    this.render();
  }

  render() {
    $(this.root).empty();
    const title = ui.h1(this.title, {style: {'align-self': 'center'}});
    const grid = this.viewerGrid!.root;
    title.style.alignContent = 'center';
    grid.style.width = 'auto';
    this.root.appendChild(ui.divV([title, grid]));
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
    let splitColumns: string[][] = Array.from({length: modeMonomerCount}, (_) => []);
    modeMonomerCount--; // minus N-terminal
    for (let i = 0; i < colLength; i++) {
      currentSplitPeptide = splitPeptidesArray[i];
      nTerminal = currentSplitPeptide.pop()!; // it is guaranteed that there will be at least one element
      currentLength = currentSplitPeptide.length;
      if (currentLength !== modeMonomerCount)
        invalidIndexes.push(i);

      for (let j = 0; j < modeMonomerCount; j++)
        splitColumns[j].push(j < currentLength ? currentSplitPeptide[j] : '-');

      splitColumns[modeMonomerCount].push(nTerminal);
    }
    modeMonomerCount--; // minus C-terminal

    //create column names list
    const columnNames = Array.from({length: modeMonomerCount}, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1}`);
    columnNames.splice(0, 0, 'N-terminal');
    columnNames.push('C-terminal');

    // filter out the columns with the same values
    if (filter) {
      splitColumns = splitColumns.filter((positionArray, index) => {
        const isRetained = new Set(positionArray).size > 1;
        if (!isRetained)
          columnNames.splice(index, 1);

        return isRetained;
      });
    }

    return splitPeptidesArray;
  }
}
