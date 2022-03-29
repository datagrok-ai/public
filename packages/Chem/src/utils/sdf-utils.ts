import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

export function _saveAsSdf() {
  //todo: load OpenChemLib (or use RDKit?)
  //todo: open dialog
  //todo: UI for choosing structure column if necessary
  //todo: UI for choosing columns with properties

  const table = grok.shell.t;
  const structureColumn = table.columns.bySemType('Molecule');
  if (structureColumn == null)
    return;

  let result = '';

  for (let i = 0; i < table.rowCount; i++) {
    try {
      const molecule: string = structureColumn.get(i);
      const mol = molecule.includes('M  END') ? OCL.Molecule.fromMolfile(molecule) : OCL.Molecule.fromSmiles(molecule);
      result += i == 0 ? '' : '\n';
      result += `${mol.toMolfile()}\n`;

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
