/* eslint-disable max-len */
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
    const notationInput = ui.input.choice('Notation', {value: null, nullable: true,
      items: [DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock], tooltipText: `Molecule notation for selected molecular column. If not specified, the exporter will try to keep existing MolBlocks as is, and convert other formats to MolBlock V2000`});
    sdfDialog.add(notationInput);
    const selectedColsInput = ui.input.bool('Selected Columns Only', {value: false, tooltipText: 'Include only the currently selected columns.'});
    sdfDialog.add(selectedColsInput);
    const visibleColsInput = ui.input.bool('Visible Columns Only', {value: true, tooltipText: 'Include only columns that are visible in the view.'});
    sdfDialog.add(visibleColsInput);
    const filteredRowsInput = ui.input.bool('Filtered Rows Only', {value: false, tooltipText: 'Include only rows that match the current filter. Can be combined with [selectedRowsOnly].'});
    sdfDialog.add(filteredRowsInput);
    const selectedRowsInput = ui.input.bool('Selected Rows Only', {value: false, tooltipText: 'Include only the currently selected rows. Can be combined with [filteredRowsOnly].'});
    sdfDialog.add(selectedRowsInput);
    if ('initFromLocalStorage' in sdfDialog && typeof sdfDialog['initFromLocalStorage'] === 'function') // temporary compatibility fix for platform versions <1.26
      sdfDialog.initFromLocalStorage();
    sdfDialog.initDefaultHistory();
    colsInput.value = cols[0];
    sdfDialog.onOK(async () => {
      const progressIndicator = DG.TaskBarProgressIndicator.create('Saving as SDF...');
      const structureColumn = colsInput.value;
      const exportOptions: DG.CsvExportOptions & {molBlockFormat?: DG.chem.Notation} = {
        molBlockFormat: notationInput.value ?? undefined,
        selectedColumnsOnly: selectedColsInput.value,
        visibleColumnsOnly: visibleColsInput.value,
        filteredRowsOnly: filteredRowsInput.value,
        selectedRowsOnly: selectedRowsInput.value,
      };
      await _saveAsSdf(table, structureColumn!, exportOptions);
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

export async function getSdfString(
  table: DG.DataFrame,
  structureColumn: DG.Column, // non-null
  exportOptions: DG.CsvExportOptions & {molBlockFormat?: DG.chem.Notation} = {selectedColumnsOnly: false, visibleColumnsOnly: true, filteredRowsOnly: false, selectedRowsOnly: false},
): Promise<string> {
  let result = '';
  //console.time('SDF export');
  // const rowIndexes = Array.from(table.rows.indexes({onlyFiltered: exportOptions.filteredRowsOnly, onlySelected: exportOptions.selectedRowsOnly}));
  // table.rows.
  const mask = DG.BitSet.create(table.rowCount);
  mask.setAll(true);
  if (exportOptions.filteredRowsOnly)
    mask.and(table.filter);
  if (exportOptions.selectedRowsOnly)
    mask.and(table.selection);

  const maskedTable = mask.anyFalse ? table.clone(mask) : table;

  const grid = grok.shell.tv != null && grok.shell.tv.dataFrame === table ? grok.shell.tv.grid :
    table.name && grok.shell.getTableView(table.name) ? grok.shell.getTableView(table.name).grid : null;
  const selectedCols = Array.from(exportOptions.selectedColumnsOnly ? table.columns.selected : table.columns);
  const visibleCols = (exportOptions.visibleColumnsOnly && grid != null ? selectedCols.filter((col) => {
    const gridCol = grid.columns.byName(col.name);
    if (gridCol == null)
      return false;
    return gridCol.settings && gridCol.settings['isPinned'] === true ?
      true : gridCol.visible; // this is needed because of how pinned columns are implemented now
  }) : selectedCols).map((col) => maskedTable.col(col.name)!);

  // precompute the molecules conversions using RDKIT service workers (much faster)

  const rdkitService = await chemCommonRdKit.getRdKitService();
  const maskedStructureCol = maskedTable.getCol(structureColumn.name);
  if (maskedStructureCol == null)
    throw new Error(`Column ${structureColumn.name} not found`);

  // if exportOptions.molBlockFormat is not set, do not convert molblocks if all of them are already in any of molblock formats
  let molBlockFormat = exportOptions.molBlockFormat;
  let dontTouchMolblocks = false;
  if (molBlockFormat == null) {
    dontTouchMolblocks = true;
    molBlockFormat = DG.chem.Notation.MolBlock;
  }

  //const molBlockFormat = exportOptions?.molBlockFormat ?? DG.chem.Notation.MolBlock;
  const molList = maskedStructureCol.toList();
  const molBlockCheck = dontTouchMolblocks ? ((mol: string) => DG.chem.isMolBlock(mol)) : (molBlockFormat === DG.chem.Notation.MolBlock ? (mol: string) => DG.chem.isMolBlock(mol) && mol.includes('V2000') :
    (mol: string) => DG.chem.isMolBlock(mol) && mol.includes('V3000'));

  const skipConversion = molList.every((mol) => !mol?.trim() || molBlockCheck(mol));
  const convertedStructures = skipConversion ? molList : await rdkitService.convertMolNotation(molList, molBlockFormat);
  const otherMolCols: {[key: string]: string[]} = {};
  for (const col of visibleCols) {
    if (col !== maskedStructureCol && col.semType === DG.SEMTYPE.MOLECULE)
      otherMolCols[col.name] = await rdkitService.convertMolNotation(col.toList(), DG.chem.Notation.Smiles);
  }

  const rowCount = maskedTable.rowCount;
  for (let rowIdx = 0; rowIdx < rowCount; rowIdx++) {
    // const molecule: string = structureColumn.get(rowIdx);
    // const mol = DG.chem.isMolBlock(molecule) ? molecule :
    //   _convertMolNotation(molecule, DG.chem.Notation.Unknown, exportOptions?.molBlockFormat ?? DG.chem.Notation.MolBlock, chemCommonRdKit.getRdKitModule());
    const mol = maskedStructureCol.isNone(rowIdx) ? '' : convertedStructures[rowIdx] ?? '';
    result += `${mol}\n`;

    // properties
    for (const col of visibleCols) {
      if (col !== maskedStructureCol) {
        let cellValue: string = '';
        if (col.semType === DG.SEMTYPE.MOLECULE)
          cellValue = otherMolCols[col.name][rowIdx] ?? '';
        else
          cellValue = col.isNone(rowIdx) ? '' : [DG.COLUMN_TYPE.QNUM, DG.COLUMN_TYPE.FLOAT].includes(col.type as DG.COLUMN_TYPE) ? DG.format(col.get(rowIdx), col.meta.format ?? undefined) : col.get(rowIdx);
        result += `>  <${col.name?.toLowerCase()?.trim() === 'molecule' ? table.columns.getUnusedName(`${col.name} (1)`) : col.name}>\n${cellValue}\n\n`;
      }
    }
    result += '$$$$\n';
  }
  //console.timeEnd('SDF export');
  return result;
}

export async function _saveAsSdf(
  table: DG.DataFrame,
  structureColumn: DG.Column,
  exportOptions: DG.CsvExportOptions & {molBlockFormat?: DG.chem.Notation} = {selectedColumnsOnly: false, visibleColumnsOnly: true, filteredRowsOnly: false, selectedRowsOnly: false},
): Promise<void> {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: UI for choosing columns with properties

  if (structureColumn == null)
    return;

  const result = await getSdfString(table, structureColumn, exportOptions);
  DG.Utils.download(table.name + '.sdf', result);
}
