import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as OCL from 'openchemlib/full';
import $ from 'cash-dom';

// import {Subscription} from 'rxjs';

// let dlg: DG.Dialog | null = null;
// let dialogSubs: Subscription[] = [];


/**  Dialog for SDF file exporter */
export function saveAsSdfDialog() {
  // here one should propose a choice of semtype, or pass it as a parameter from
  // ui
  // const structureColumn = table.columns.bySemType('Molecule');
  const table = grok.shell.t;
  const moleculeCols = table.columns.bySemTypeAll('Molecule');
  const macromoleculeCols = table.columns.bySemTypeAll('Macromolecule');
  if (moleculeCols.length === 0 && macromoleculeCols.length === 0) {
    const dlg = ui.dialog({title: 'Save as SDF'})
      .add(
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
    $(dlg.getButton('CANCEL')).hide();
    dlg.show({x: 350, y: 100});
  } else {
    const cols = moleculeCols.concat(macromoleculeCols);
    const colsInput = ui.choiceInput('Choose column:', cols[0], cols);
    const dlg = ui.dialog({title: 'Save as SDF'})
      .add(ui.div([
        colsInput.root,
      ]))
      .onOK(() => {
        const structureColumn = colsInput.value;
        _saveAsSdf(table, structureColumn!);
      });
    dlg.show({x: 350, y: 100});

    // dialogSubs.push(dlg.onClose.subscribe((value) => {
    //   dialogSubs.forEach((s) => {s.unsubscribe();});
    //   dialogSubs = [];
    }));
  }
}

export function _saveAsSdf(
  table: DG.DataFrame,
  structureColumn: DG.Column,
): void {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing structure column if necessary
  //todo: UI for choosing columns with properties

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
}
