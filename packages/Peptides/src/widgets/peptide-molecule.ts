import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';

/**
 * 3D representation widget of peptide molecule.
 *
 * @export
 * @param {string} pep Peptide string.
 * @return {Promise<DG.Widget>} Widget.
 */
export async function peptideMoleculeWidget(pep: string, currentTable: DG.DataFrame): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating NGL view');
  const separator = currentTable.columns.bySemType(C.SEM_TYPES.MACROMOLECULE)!.tags[C.TAGS.SEPARATOR];

  let widgetHost;
  let smiles = '';
  let molfileStr = '';
  try {
    try {
      const params = {table: currentTable};
      const result = await grok.functions.call('Customerextensions:getPeptideStructure', params) as string[];
      if (result.length !== 0) {
        smiles = result[0];
        molfileStr = result[1];
        throw new Error(`Found structure in DB`);
      }

      smiles = getMolecule(pep, separator);
      if (smiles == '')
        throw new Error('Couldn\'t get smiles');

      molfileStr = (await grok.functions.call('Peptides:SmiTo3D', {smiles})) as string;
    } catch (e) {
      console.warn(e);
    }

    try {
      molfileStr = molfileStr.replaceAll('\\n', '\n');
      const stringBlob = new Blob([molfileStr], {type: 'text/plain'});
      const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});

      //@ts-ignore
      const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
      //@ts-ignore
      stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
        stage.setSize(300, 300);
        comp.addRepresentation('ball+stick');
        comp.autoView();
      });
      const sketch = grok.chem.svgMol(molfileStr);
      const panel = ui.divH([sketch]);

      widgetHost = ui.div([panel, nglHost]);
    } catch (e) {
      widgetHost = ui.divText('Couldn\'t get peptide structure');
    }
  } catch (e) {
    widgetHost = ui.divText('Couldn\'t get peptide structure');
  }
  pi.close();
  return new DG.Widget(widgetHost);
}

//FIXME: doesn't work after removing chemPalette
export function getMolecule(pep: string, separator: string): string {
  const split = pep.split(separator);
  const mols = [];
  //@ts-ignore
  const chemPalette = PeptidesModel.chemPalette;
  for (let i = 1; i < split.length - 1; i++) {
    if (split[i] in chemPalette.AASmiles) {
      const aar = chemPalette.AASmiles[split[i]];
      mols[i] = aar.substring(0, aar.length - 1);
    } else if (!split[i] || split[i] == '-')
      mols[i] = '';
    else
      return '';
  }

  return mols.join('') + 'O';
}
