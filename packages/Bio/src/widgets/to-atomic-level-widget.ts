import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package, PackageFunctions} from '../package';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';


export async function toAtomicLevelSingle(sequence: DG.SemanticValue,
  range?: [number, number]): Promise<{mol: string, errorText: string, scrollableLength?: number}> {
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
    const seqHelper = await PackageFunctions.getSeqHelper();
    const seqSh = seqHelper.getSeqHandler(sequence.cell.column);
    if (!seqSh) {
      errorText = 'No sequence handler found';
      return {errorText, mol: ''};
    }

    const seqSplitted = seqSh.getSplitted(sequence.cell.rowIndex);
    let maxLength = 35;
    if (seqSh.isHelm() && !seqSplitted.graphInfo?.polymerTypes?.some((pt) => pt !== 'RNA'))
      maxLength = 105;

    // helm and biln will have cyclization marks, so we need to use POM to convert them
    const shouldUsePOM = (seqSplitted.graphInfo?.connections?.length ?? 0) > 0 || seqSh.units === NOTATION.CUSTOM;
    const isHelmWithMultiplePolymerTypes = seqSh.isHelm() &&
      (new Set((seqSplitted.graphInfo?.polymerTypes ?? []))).size > 1;

    if (!range && (seqSplitted?.length ?? 100) > maxLength) {
      errorText = 'Maximum number of monomers for molecular conversion is ' + maxLength;
      // fasta/separator sequences are single linear chains, so they can still be viewed in scrollable windows;
      // isSeparator() is true for BILN columns too (they carry a separator tag), so check the notation strictly
      // and exclude columns with a custom notation provider
      const canScroll = [NOTATION.FASTA, NOTATION.SEPARATOR].includes(seqSh.notation) &&
        !sequence.cell.column.temp[SeqTemps.notationProvider];
      return canScroll ? {errorText, mol: '', scrollableLength: seqSplitted.length} : {errorText, mol: ''};
    }
    let singleValCol = DG.Column.fromStrings('singleVal', [sequence.value]);
    // copy over all the tags
    Object.entries(sequence.cell.column.tags).forEach(([key, value]) => {
      singleValCol.setTag(key, value as string);
    });

    // if column has notation provider, we need to copy it over
    const notationProvider = sequence.cell.column.temp[SeqTemps.notationProvider];
    if (notationProvider)
      singleValCol.temp[SeqTemps.notationProvider] = notationProvider;

    if (range) {
      singleValCol = seqHelper.getSeqHandler(singleValCol).getRegion(range[0], range[1], 'singleVal');
      if (notationProvider)
        singleValCol.temp[SeqTemps.notationProvider] = notationProvider;
    }
    const sDf = DG.DataFrame.fromColumns([singleValCol]);

    await PackageFunctions.toAtomicLevel(sDf, singleValCol,
      shouldUsePOM || isHelmWithMultiplePolymerTypes, false);
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

function drawMoleculeHost(mol: string, sequence: DG.SemanticValue): HTMLElement {
  let width = 300;
  let height = 300;
  const tagW = Number.parseInt(sequence.cell.column.getTag('.toAtomicWidgetWidth') ?? '');
  const tagH = Number.parseInt(sequence.cell.column.getTag('.toAtomicWidgetHeight') ?? '');
  if (tagW && Number.isFinite(tagW))
    width = tagW;
  if (tagH && Number.isFinite(tagH))
    height = tagH;

  const root = grok.chem.drawMolecule(mol, width, height, false);
  root.style.cursor = 'pointer';
  ui.tooltip.bind(root, 'Click to expand');
  root.onclick = () => {
    const width = window.innerWidth - 200;
    const height = window.innerHeight - 200;
    const bigMol = grok.chem.drawMolecule(mol, width, height, false);
    ui.dialog({title: 'Molecule'}).add(bigMol).showModal(true);
  };
  return root;
}

const SCROLL_WINDOW_SIZE = 35;
const SCROLL_RENDER_DEBOUNCE_MS = 300;

/** For long linear single-chain sequences: shows the structure of a {@link SCROLL_WINDOW_SIZE}-monomer
 * window with a slider to scroll along the sequence. Positions in the UI are 1-based. */
function toAtomicLevelScrollableWidget(sequence: DG.SemanticValue, seqLength: number): DG.Widget {
  const maxStart = seqLength - SCROLL_WINDOW_SIZE + 1; // 1-based
  const startLabel = ui.divText('');
  const endLabel = ui.divText('');
  const molHost = ui.div([ui.loader()]);
  let renderId = 0;
  let renderTimer: number | null = null;
  const offsetEl = ui.span([]);
  const updateLabels = (start: number) => {
    startLabel.textContent = `${start}`;
    endLabel.textContent = `${Math.min(start + SCROLL_WINDOW_SIZE - 1, seqLength)}`;
    offsetEl.textContent = `${start}`;
  };

  const render = async (start: number) => {
    const currentId = ++renderId;
    const res = await toAtomicLevelSingle(sequence,
      [start - 1, Math.min(start + SCROLL_WINDOW_SIZE - 1, seqLength) - 1]);
    if (currentId !== renderId)
      return;
    ui.empty(molHost);
    molHost.append(res.errorText || !res.mol ?
      ui.divText(res.errorText || 'No structure generated') : drawMoleculeHost(res.mol, sequence));
  };
  const slider = ui.input.slider('', {value: 1, min: 1, max: maxStart, step: 1,
    tooltipText: 'Sequence offset to translate',
    showPlusMinus: true,
    onValueChanged: (value) => {
      const start = Math.min(Math.max(Math.round(value ?? 1), 1), maxStart);
      updateLabels(start);
      if (renderTimer !== null)
        window.clearTimeout(renderTimer);
      renderTimer = window.setTimeout(() => {
        renderTimer = null;
        render(start);
      }, SCROLL_RENDER_DEBOUNCE_MS);
    }});
  slider.input.style.width = '100%';
  slider.addOptions(offsetEl);
  updateLabels(1);
  render(1);
  return DG.Widget.fromRoot(ui.divV([
    ui.divH([startLabel, endLabel], {style: {justifyContent: 'space-between'}}),
    molHost,
    slider.root,
  ]));
}

export async function toAtomicLevelWidget(sequence: DG.SemanticValue): Promise<DG.Widget> {
  const res = await toAtomicLevelSingle(sequence);
  if (res.scrollableLength)
    return toAtomicLevelScrollableWidget(sequence, res.scrollableLength);
  if (res.errorText || !res.mol)
    return DG.Widget.fromRoot(ui.divText(res.errorText ?? 'No structure generated'));
  try {
    const molSemanticValue = DG.SemanticValue.fromValueType(res.mol, DG.SEMTYPE.MOLECULE);
    const panel = ui.panels.infoPanel(molSemanticValue);
    let molPanel: DG.Widget | null = null;
    if (panel) {
      const acc = ui.accordion('Sequence Molfile details');
      acc.addPane('Explore', () => panel.root);
      molPanel = DG.Widget.fromRoot(acc.root);
    }

    const root = drawMoleculeHost(res.mol, sequence);
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
      molBlock3D = (await grok.functions.call('Bio:Embed', {molecule: molBlock2D})) as unknown as string;
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
