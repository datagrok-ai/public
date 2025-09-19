import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package, PackageFunctions} from '../package';
import * as api from '../package-api';

export async function toAtomicLevelSingle(sequence: DG.SemanticValue): Promise<{mol: string, errorText: string}> {
  let errorText = '';
  try {
    if (!sequence || !sequence.value) {
      errorText = 'No sequence provided';
      return {errorText, mol: ''};
    }
    if (!sequence.cell || !sequence.cell.dart || !sequence.cell.dataFrame || !sequence.cell.column) {
      errorText = 'Atomic level conversion requeires a sequence column';
      return {errorText, mol: ''};
    }
    const supportedUnits: string[] = [NOTATION.FASTA, NOTATION.SEPARATOR, NOTATION.HELM];
    //todo: add support for custom notations
    if (!supportedUnits.includes(sequence.cell.column.meta.units?.toLowerCase() ?? '')) {
      errorText = 'Unsupported sequence notation. please use Bio | Polytool | Convert';
      return {errorText, mol: ''};
    }
    const seqHelper = await PackageFunctions.getSeqHelper();
    const seqSh = seqHelper.getSeqHandler(sequence.cell.column);
    if (!seqSh) {
      errorText = 'No sequence handler found';
      return {errorText, mol: ''};
    }
    if ((seqSh.getSplitted(sequence.cell.rowIndex, 50)?.length ?? 100) > 40) {
      errorText = 'Maximum number of monomers is 40';
      return {errorText, mol: ''};
    }
    const singleValCol = DG.Column.fromStrings('singleVal', [sequence.value]);
    const sDf = DG.DataFrame.fromColumns([singleValCol]);
    // copy over all the tags
    Object.entries(sequence.cell.column.tags).forEach(([key, value]) => {
      singleValCol.setTag(key, value as string);
    });
    await PackageFunctions.toAtomicLevel(sDf, singleValCol, sequence.cell.column.meta.units === NOTATION.HELM, false);
    if (sDf.columns.length < 2) {
      errorText = 'No structure generated';
      return {errorText, mol: ''};
    }
    const molCol = sDf.columns.byIndex(1);
    const molfile = molCol.get(0);
    if (!molfile) {
      errorText = 'No structure generated';
      return {errorText, mol: ''};
    }
    return {errorText: '', mol: molfile as string};
  } catch (e) {
    _package.logger.error(e);
  }

  errorText = 'No Structure generated';
  return {errorText, mol: ''};
}

export async function toAtomicLevelWidget(sequence: DG.SemanticValue): Promise<DG.Widget> {
  const res = await toAtomicLevelSingle(sequence);
  if (res.errorText || !res.mol)
    return DG.Widget.fromRoot(ui.divText(res.errorText ?? 'No structure generated'));
  try {
    const molSemanticValue = DG.SemanticValue.fromValueType(res.mol, DG.SEMTYPE.MOLECULE);
    const panel = ui.panels.infoPanel(molSemanticValue);
    let molPanel: DG.Widget | null = null;
    if (panel)
      molPanel = DG.Widget.fromRoot(panel.root);


    const root = grok.chem.drawMolecule(res.mol, 300, 300, false);
    root.style.cursor = 'pointer';
    ui.tooltip.bind(root, 'Click to expand');
    root.onclick = () => {
      const width = window.innerWidth - 200;
      const height = window.innerHeight - 200;
      const bigMol = grok.chem.drawMolecule(res.mol, width, height, false);
      ui.dialog({title: 'Molecule'}).add(bigMol).showModal(true);
    };
    if (molPanel)
      molPanel.root.prepend(root);
    return molPanel ?? DG.Widget.fromRoot(root);
  } catch (e) {
    _package.logger.error(e);
  }
  return DG.Widget.fromRoot(ui.divText('No structure generated'));
}

/**
 * 3D representation widget of macromolecule.
 *
 * @export
 * @return {Promise<DG.Widget>} Widget.
 */
export async function molecular3DStructureWidget(
  sequence: DG.SemanticValue
): Promise<DG.Widget> {
  const pi = DG.TaskBarProgressIndicator.create('Creating 3D view');
  let widgetHost;
  let molBlock3D = '';
  try {
    // make sure biostructure viewer package is loaded.
    await DG.Func.find({name: 'getPdbHelper'})[0]?.apply({});
    try {
      const result = await toAtomicLevelSingle(sequence);//await getMacroMol(atomicCodes!);
      if (result.errorText || !result.mol) {
        widgetHost = ui.divText(result.errorText ?? 'No structure generated');
        pi.close();
        return new DG.Widget(widgetHost);
      }
      const molBlock2D = result.mol;
      molBlock3D = await api.scripts.embed(molBlock2D)as unknown as string;
      // rdfkit sometimes fails to convert molv3 to molv2, so we try to convert it via the OCL
      const OCLMol = OCL.Molecule.fromMolfile(molBlock3D);
      if (OCLMol)
        molBlock3D = OCLMol.toMolfile();
      else
        console.warn('Failed to convert molv3 to molv2');

      //molBlock3D = grok.chem.convert(molBlock3D, grok.chem.Notation.Unknown, grok.chem.Notation.MolBlock);
    } catch (e) {
      console.warn(e);
    }

    try {
      molBlock3D = molBlock3D.replaceAll('\\n', '\n');
      const stringBlob = new Blob([molBlock3D], {type: 'text/plain'});
      const nglHost = ui.div([], {classes: 'd4-ngl-viewer', id: 'ngl-3d-host'});
      nglHost.style.setProperty('height', '100%', 'important');
      //@ts-ignore
      const stage = new NGL.Stage(nglHost, {backgroundColor: 'white'});
      //@ts-ignore
      stage.loadFile(stringBlob, {ext: 'sdf'}).then(function(comp: NGL.StructureComponent) {
        stage.setSize(300, 300);
        comp.addRepresentation('ball+stick');
        comp.autoView();
      });
      widgetHost = ui.div([nglHost], {style: {aspectRatio: '1'}});
    } catch (e) {
      widgetHost = ui.divText('Couldn\'t get 3D structure');
    }
  } catch (e) {
    widgetHost = ui.divText('Couldn\'t get 3D structure');
  }
  pi.close();
  return new DG.Widget(widgetHost);
}
