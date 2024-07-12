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
  else if (cols.length === 1)
    _saveAsSdf(table, cols[0]);
  else {
    const sdfDialog = ui.dialog({title: 'Save as SDF'});
    sdfDialog.root.style.width = '250px';
    const colsChoiceDF = DG.DataFrame.fromColumns(cols);
    const colsInput = ui.input.column('Molecules', {table: colsChoiceDF, value: cols[0]});

    sdfDialog.add(colsInput)
      .onOK(() => {
        const structureColumn = colsInput.value;
        _saveAsSdf(table, structureColumn!);
      });
    if (cols.length > 1) {
      const text = ui.divText(`Other ${DG.SEMTYPE.MOLECULE} colums saved as SMILES`);
      sdfDialog.add(text);
      $(text).css('color', 'gray');
      $(text).css('font-size', '12px');
    }
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
): string {
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    const molecule: string = structureColumn.get(i);
    const mol = DG.chem.isMolBlock(molecule) ? molecule :
      _convertMolNotation(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.MolBlock,
        chemCommonRdKit.getRdKitModule());
    result += `${mol}\n`;

    // properties
    for (const col of table.columns) {
      if (col !== structureColumn) {
        let cellValue = col.get(i);
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
): void {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: UI for choosing columns with properties

  if (structureColumn == null)
    return;

  const result = getSdfString(table, structureColumn);
  DG.Utils.download(table.name + '.sdf', result);
}
