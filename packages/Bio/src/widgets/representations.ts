import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getMolfilesFromSingleSeq, HELM_CORE_LIB_FILENAME} from '../utils/utils';
import {getMacroMol} from '../utils/atomic-works';

/**
 * 3D representation widget of macromolecule.
 *
 * @export
 * @param {DG.Cell} macroMolecule macromolecule cell.
 * @return {Promise<DG.Widget>} Widget.
 */
export async function representationsWidget(macroMolecule: DG.Cell, monomersLibObject: any[]): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

  let widgetHost;
  let molBlock3D = '';
  try {
    try {
      const atomicCodes = getMolfilesFromSingleSeq(macroMolecule, monomersLibObject);
      const result = await getMacroMol(atomicCodes!);
      const molBlock2D = result[0];
      molBlock3D = (await grok.functions.call('Bio:Embed', {molBlock2D})) as string;
    } catch (e) {
      console.warn(e);
    }

    try {
      molBlock3D = molBlock3D.replaceAll('\\n', '\n');
      const stringBlob = new Blob([molBlock3D], {type: 'text/plain'});
      const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});

      //@ts-ignore
      const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
      //@ts-ignore
      stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
        stage.setSize(300, 300);
        comp.addRepresentation('ball+stick');
        comp.autoView();
      });
      const sketch = grok.chem.svgMol(molBlock3D);
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
