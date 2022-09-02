import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as OCL from 'openchemlib/full';
import $ from 'cash-dom';

/**  Dialog for SDF file exporter */
export function saveAsSdfDialog() {
  const table = grok.shell.t;
  let cols = table.columns.bySemTypeAll('Molecule');
  cols = cols.concat(table.columns.bySemTypeAll('Macromolecule'));
  if (cols.length === 0)
    grok.shell.warning('This table does not contain Molecule columns, unable to save as SDF');
  else if (cols.length === 1)
    _saveAsSdf(table, cols[0]);
  else {
    const sdfDialog = ui.dialog({title: 'Save as SDF'});
    sdfDialog.root.style.width = '250px';
    const colsChoiceDF = DG.DataFrame.fromColumns(cols);
    const colsInput = ui.columnInput('Molecules', colsChoiceDF, cols[0]);

    sdfDialog.add(colsInput)
      .onOK(() => {
        const structureColumn = colsInput.value;
        _saveAsSdf(table, structureColumn!);
      });
    sdfDialog.show({x: 350, y: 100});
    $(sdfDialog.root).find('.grok-font-icon-help').remove();
    $(sdfDialog.root).find('.ui-input-label').css('max-width', 'fit-content');
  }
}

export function _saveAsSdf(
  table: DG.DataFrame,
  structureColumn: DG.Column,
): void {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing columns with properties

  const pi = DG.TaskBarProgressIndicator.create('Saving as SDF...');

  if (structureColumn == null)
    return;

  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      const molecule: string = structureColumn.get(i);
      const mol = molecule.includes('M  END') ? molecule : OCL.Molecule.fromSmiles(molecule).toMolfile();
      result += i == 0 ? '' : '\n';
      result += `${mol}\n`;

      // properties
      for (const col of table.columns) {
        if (col !== structureColumn)
          result += `>  <${col.name}>\n${col.get(i)}\n\n`;
      }

      result += '$$$$';
    } catch (error) {
      console.error(error);
    }
  }

  DG.Utils.download(table.name + '.sdf', result);
  pi.close();
}
