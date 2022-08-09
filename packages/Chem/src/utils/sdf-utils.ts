import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as OCL from 'openchemlib/full';
import $ from 'cash-dom';

import {Subscription} from 'rxjs';

let sdfDialog: DG.Dialog | null = null;
let sdfDialogSubs: Subscription[] = [];

/**  Dialog for SDF file exporter */
export function saveAsSdfDialog() {
  const table = grok.shell.t;
  let cols = table.columns.bySemTypeAll('Molecule');
  cols = cols.concat(table.columns.bySemTypeAll('Macromolecule'));
  if (cols.length === 1)
    _saveAsSdf(table, cols[0]);
  else if (sdfDialog === null) {
    sdfDialog = ui.dialog({title: 'Save as SDF'});
    if (cols.length === 0) {
      sdfDialog.add(
        ui.divText(
          'This table does not contain Molecule/Macromolecule columns, unable to save as SDF',
          {style:
            {
              'box-sizing': 'border-box',
              'width': '200px',
              'padding': '5px',
              'display': 'inline-block',
              'text-align': 'center',
            },
          },
        ),
      )
        .onOK(() => {});
      $(sdfDialog.getButton('CANCEL')).hide();
    } else {
      const colsInput = ui.choiceInput('Choose column:', cols[0], cols);
      sdfDialog.add(ui.div([
        colsInput.root,
      ]))
        .onOK(() => {
          const structureColumn = colsInput.value;
          _saveAsSdf(table, structureColumn!);
        });
    }
    sdfDialog.show({x: 350, y: 100});
    sdfDialogSubs.push(sdfDialog.onClose.subscribe((value) => {
      sdfDialogSubs.forEach((s) => {s.unsubscribe();});
      sdfDialogSubs = [];
      sdfDialog = null;
    }));
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
