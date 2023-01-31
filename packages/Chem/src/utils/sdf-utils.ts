import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as chemCommonRdKit from './chem-common-rdkit';
import {MolNotation, _convertMolNotation} from './convert-notation-utils';
import {isMolBlock} from './convert-notation-utils';
import $ from 'cash-dom';

/**  Dialog for SDF file exporter */
export async function saveAsSdfDialog() {
  const table = grok.shell.t;
  const cols = table.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
  if (cols.length === 0)
    grok.shell.warning(`This table does not contain ${DG.SEMTYPE.MOLECULE} columns, unable to save as SDF`);
  else if (cols.length === 1)
    await _saveAsSdf(table, cols[0]);
  else {
    const sdfDialog = ui.dialog({title: 'Save as SDF'});
    sdfDialog.root.style.width = '250px';
    const colsChoiceDF = DG.DataFrame.fromColumns(cols);
    const colsInput = ui.columnInput('Molecules', colsChoiceDF, cols[0]);

    sdfDialog.add(colsInput)
      .onOK(async () => {
        const structureColumn = colsInput.value;
        await _saveAsSdf(table, structureColumn!);
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

export async function getSdfString(
  table: DG.DataFrame,
  structureColumn: DG.Column, // non-null
): Promise<string> {
  let result = '';
  for (let i = 0; i < table.rowCount; i++) {
    const molecule: string = structureColumn.get(i);
    const mol = isMolBlock(molecule) ? molecule :
      _convertMolNotation(molecule, MolNotation.Unknown, MolNotation.MolBlock, chemCommonRdKit.getRdKitModule());
    result += i == 0 ? '' : '\n';
    result += `${mol}\n`;

    // properties
    for (const col of table.columns) {
      if (col !== structureColumn) {
        let cellValue = col.get(i);
        // convert to SMILES if necessary
        if (col.semType === DG.SEMTYPE.MOLECULE) {
          cellValue = _convertMolNotation(cellValue, MolNotation.Unknown, 
            MolNotation.Smiles, chemCommonRdKit.getRdKitModule());
        }
        result += `>  <${col.name}>\n${cellValue}\n\n`;
      }
    }
    result += '$$$$';
  }
  return result;
}

export async function _saveAsSdf(
  table: DG.DataFrame,
  structureColumn: DG.Column,
): Promise<void> {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: UI for choosing columns with properties

  if (structureColumn == null)
    return;

  const result = await getSdfString(table, structureColumn);
  DG.Utils.download(table.name + '.sdf', result);
}
