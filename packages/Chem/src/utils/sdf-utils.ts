import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as chemCommonRdKit from './chem-common-rdkit';
import {_convertMolNotation, convertNotationForColumn} from './convert-notation-utils';
import $ from 'cash-dom';

/**  Dialog for SDF file exporter */
export function saveAsSdfDialog() {
  const table = grok.shell.t;
  const cols = table.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
  if (cols.length === 0)
    grok.shell.warning(`This table does not contain ${DG.SEMTYPE.MOLECULE} columns, unable to save as SDF`);
  else {
    const sdfDialog = ui.dialog({title: 'Save as SDF'});
    const colsChoiceDF = DG.DataFrame.fromColumns(cols);
    const colsInput = ui.input.column('Molecules', {table: colsChoiceDF, value: cols[0], nullable: false,
      filter: (col) => col.semType === DG.SEMTYPE.MOLECULE});
    sdfDialog.add(colsInput);
    if (cols.length > 1) {
      const text = ui.divText(`Other ${DG.SEMTYPE.MOLECULE} colums saved as SMILES`);
      sdfDialog.add(text);
      $(text).css('color', 'gray');
      $(text).css('font-size', '12px');
    }
    const selectedColsInput = ui.input.bool('Selected Columns Only', {value: false});
    sdfDialog.add(selectedColsInput);
    const visibleColsInput = ui.input.bool('Visible Columns Only', {value: true});
    sdfDialog.add(visibleColsInput);
    const filteredRowsInput = ui.input.bool('Filtered Rows Only', {value: false});
    sdfDialog.add(filteredRowsInput);
    const selectedRowsInput = ui.input.bool('Selected Rows Only', {value: false});
    sdfDialog.add(selectedRowsInput);
    sdfDialog.initFromLocalStorage();
    colsInput.value = cols[0];
    sdfDialog.onOK(() => {
      const progressIndicator = DG.TaskBarProgressIndicator.create('Saving as SDF...');
      const structureColumn = colsInput.value;
      const exportOptions: DG.CsvExportOptions = {
        selectedColumnsOnly: selectedColsInput.value,
        visibleColumnsOnly: visibleColsInput.value,
        filteredRowsOnly: filteredRowsInput.value,
        selectedRowsOnly: selectedRowsInput.value,
      };
      _saveAsSdf(table, structureColumn!, exportOptions);
      progressIndicator.close();
    });
    sdfDialog.show({x: 350, y: 100});
    $(sdfDialog.root).find('.grok-font-icon-help').remove();
    $(sdfDialog.root).find('.ui-input-label').css('max-width', 'fit-content');
  }
}

export async function getSdfStringAsync(structureColumn: DG.Column): Promise<string> {
  const table: DG.DataFrame = structureColumn.dataFrame;
  const convertedStruct = await convertNotationForColumn(structureColumn, DG.chem.Notation.MolBlock);
  const convertedOther = [];
  for (const col of table.columns) {
    if (col !== structureColumn) {
      if (col.semType === DG.SEMTYPE.MOLECULE)
        convertedOther.push(await convertNotationForColumn(col, DG.chem.Notation.Smiles));
    }
  }
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    result += `${convertedStruct[i]}\n`;
    let molColIdx = 0;
    for (const col of table.columns) {
      if (col.name === structureColumn.name)
        continue;

      if (col.semType === DG.SEMTYPE.MOLECULE) {
        result += `>  <${col.name}>\n${convertedOther[molColIdx]}\n\n`;
        ++molColIdx;
        continue;
      }
      result += `>  <${col.name}>\n${col.get(i)}\n\n`;
    }
    result += '$$$$\n';
  }
  return result;
}

export function getSdfString(
  table: DG.DataFrame,
  structureColumn: DG.Column, // non-null
  exportOptions: DG.CsvExportOptions = {selectedColumnsOnly: false, visibleColumnsOnly: true, filteredRowsOnly: false, selectedRowsOnly: false},
): string {
  let result = '';
  const rowIndexes = table.rows.indexes({onlyFiltered: exportOptions.filteredRowsOnly, onlySelected: exportOptions.selectedRowsOnly});
  for (const rowIdx of rowIndexes) {
    const molecule: string = structureColumn.get(rowIdx);
    const mol = DG.chem.isMolBlock(molecule) ? molecule :
      _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock, chemCommonRdKit.getRdKitModule());
    result += `${mol}\n`;

    const grid = grok.shell.tv != null && grok.shell.tv.dataFrame === table ? grok.shell.tv.grid :
      table.name && grok.shell.getTableView(table.name) ? grok.shell.getTableView(table.name).grid : null;
    const selectedCols = Array.from(exportOptions.selectedColumnsOnly ? table.columns.selected : table.columns);
    const visibleCols = exportOptions.visibleColumnsOnly && grid != null ? selectedCols.filter((col) => {
      const gridCol = grid.columns.byName(col.name);
      if (gridCol == null)
        return false;
      return gridCol.settings && gridCol.settings['isPinned'] === true ?
        true : gridCol.visible; // this is needed because of how pinned columns are implemented now
    }) : selectedCols;

    // properties
    for (const col of visibleCols) {
      if (col !== structureColumn) {
        let cellValue = col.isNone(rowIdx) ? '' : col.get(rowIdx);
        // convert to SMILES if necessary
        if (col.semType === DG.SEMTYPE.MOLECULE) {
          cellValue = _convertMolNotation(cellValue, DG.chem.Notation.Unknown,
            DG.chem.Notation.Smiles, chemCommonRdKit.getRdKitModule());
        }
        result += `>  <${col.name}>\n${cellValue}\n\n`;
      }
    }
    result += '$$$$\n';
  }
  return result;
}

export function _saveAsSdf(
  table: DG.DataFrame,
  structureColumn: DG.Column,
  exportOptions: DG.CsvExportOptions = {selectedColumnsOnly: false, visibleColumnsOnly: true, filteredRowsOnly: false, selectedRowsOnly: false},
): void {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: UI for choosing columns with properties

  if (structureColumn == null)
    return;

  const result = getSdfString(table, structureColumn, exportOptions);
  DG.Utils.download(table.name + '.sdf', result);
}
