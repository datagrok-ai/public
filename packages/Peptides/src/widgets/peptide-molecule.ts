import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// import {ChemPalette} from '../utils/chem-palette';
import {PeptidesController} from '../peptides';

/**
 * 3D representation widget of peptide molecule.
 *
 * @export
 * @param {string} pep Peptide string.
 * @return {Promise<DG.Widget>} Widget.
 */
export async function peptideMoleculeWidget(pep: string): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating NGL view');

  const smiles = getMolecule(pep);
  if (smiles == '')
    return new DG.Widget(ui.divH([]));


  let molfileStr = (await grok.functions.call('Peptides:SmiTo3D', {smiles}));

  molfileStr = molfileStr.replaceAll('\\n', '\n'); ;
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
  const sketch = grok.chem.svgMol(smiles);
  const panel = ui.divH([sketch]);

  pi.close();

  return new DG.Widget(ui.div([panel, nglHost]));
}

export function getMolecule(pep: string): string {
  const split = pep.split('-');
  const mols = [];
  const chemPalette = PeptidesController.chemPalette;
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
