import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getMolfilesFromSingleSeq} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {HELM_CORE_LIB_FILENAME} from '@datagrok-libraries/bio/src/utils/const';

/**
 * @export
 * @param {DG.Column} col macromolecule cell.
 * @return {Promise<DG.Widget>} Widget.
 */
export function getMacroMolColumnPropertyPanel(col: DG.Column): DG.Widget {
  const NONE = 'None';
  const scaffoldColName = 'short';

  // TODO: replace with an efficient version, bySemTypesExact won't help; GROK-8094
  const columnsList = Array.from(col.dataFrame.columns as any).filter(
    (c: any) => c.semType === DG.SEMTYPE.MOLECULE).map((c: any) => c.name);
  const columnsSet = new Set(columnsList);
  columnsSet.delete(col.name);

  const monomerWidth = ui.choiceInput('Monomer width',
    (col?.temp['monomer-width'] != null) ? col.temp['monomer-width'] : 'short',
    ['short', 'long'],
    (s: string) => {
      col.temp['monomer-width'] = s;
      col.setTag('.calculatedCellRender', '0');
      col.dataFrame.fireValuesChanged();
    });
  monomerWidth.setTooltip('In short mode, only the first character should be visible, followed by .. if there are more characters');

  const colorCode = ui.boolInput('Color code',
    (col?.temp['color-code'] != null) ? col.temp['color-code'] : true,
    (v: boolean) => {
      col.temp['color-code'] = v;
      col.dataFrame.fireValuesChanged();
    });
  colorCode.setTooltip('Color code');

  const referenceSequence = ui.stringInput('Reference sequence', (col?.temp['reference-sequence'] != null) ? col?.temp['reference-sequence'] : '', (v: string) => {
    col.temp['reference-sequence'] = v;
    col.dataFrame.fireValuesChanged();
  });
  referenceSequence.setTooltip('Reference sequence is not empty, then the sequence will be render ' + '\n' +
    'as a difference from the reference sequence');

  const compareWithCurrent = ui.boolInput('Compare with current',
    (col?.temp['compare-with-current'] != null) ? col.temp['compare-with-current'] : true,
    (v: boolean) => {
      col.temp['compare-with-current'] = v;
      col.dataFrame.fireValuesChanged();
    });
  compareWithCurrent.setTooltip('When on, all sequences get rendered in the "diff" mode');

  const rdKitInputs = ui.inputs([
    monomerWidth,
    colorCode,
    referenceSequence,
    compareWithCurrent
  ]);

  return new DG.Widget(rdKitInputs);
}

/**
 * 3D representation widget of macromolecule.
 *
 * @export
 * @param {DG.Cell} macroMolecule macromolecule cell.
 * @param {any[]} monomersLibObject
 * @return {Promise<DG.Widget>} Widget.
 */
export async function representationsWidget(macroMolecule: DG.Cell, monomersLibObject: any[]): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');

  let widgetHost;
  let molBlock3D = '';
  try {
    try {
      const atomicCodes = getMolfilesFromSingleSeq(macroMolecule, monomersLibObject);
      const result = ''//await getMacroMol(atomicCodes!);
      const molBlock2D = result[0];
      molBlock3D = (await grok.functions.call('Bio:Embed', {molBlock2D})) as unknown as string;
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
