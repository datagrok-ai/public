import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMacromoleculeAssembler} from '@datagrok-libraries/bio/src/monomer-works/assembler';
import {_package, PackageFunctions, toAtomicLevelDo} from '../package';
import {SeqTemps} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';


/** One sequence cell to a molblock. The linear layout by default: the assembler falls back per row to the
 * non-linear one for rows that are not plain chains (cycles, cross-links, CHEM monomers, mixed polymer types).
 * @param {DG.SemanticValue} sequence a cell of a Macromolecule column
 * @param {[number, number]} range a window of a long FASTA/separator sequence (0-based, inclusive)
 * @param {boolean} nonlinear the non-linear layout regardless of the row */
export async function toAtomicLevelSingle(sequence: DG.SemanticValue,
  range?: [number, number], nonlinear: boolean = false,
): Promise<{mol: string, errorText: string, scrollableLength?: number}> {
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
    // Converting and drawing 200 monomers takes about 300 ms; the full-screen dialog (a click on the molecule) shows
    // the details of a long one. An all-RNA HELM row counts sugar, base and phosphate as three positions, so 360 is
    // about 120 nucleotides
    let maxLength = 120;
    if (seqSh.isHelm() && !seqSplitted.graphInfo?.polymerTypes?.some((pt) => pt !== 'RNA'))
      maxLength = 360;

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

    // Linear unless asked otherwise: rows with connections, custom notations or mixed polymer types fall back to the
    // non-linear layout in the assembler, per row. No warning balloon: the panel shows the error itself
    await toAtomicLevelDo(sDf, singleValCol, nonlinear, false, false);
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
    showMoleculeDialog(mol).catch((err) => _package.logger.error(err));
  };
  return root;
}

/** The full-screen dialog draws the molecule at this many times the resolution it shows it at */
const DIALOG_RESOLUTION = 2;
/** At most this many canvas pixels for that drawing (4 bytes each; Chem draws it uncached, so the memory goes with
 * the dialog) */
const DIALOG_MAX_PIXELS = 25_000_000;
const LENS_MAGNIFICATION = 3;
/** The lens' width and height, as a fraction of the dialog's molecule view */
const LENS_FRACTION = 1 / 3;

/** The molecule in a full-screen dialog, drawn at {@link DIALOG_RESOLUTION} times the shown resolution, with a
 * magnifying lens under the mouse: the part around the pointer {@link LENS_MAGNIFICATION} times larger. */
export async function showMoleculeDialog(mol: string): Promise<DG.Dialog> {
  const width = Math.max(200, window.innerWidth - 200);
  const height = Math.max(100, window.innerHeight - 200);
  const view = ui.div([], {style: {position: 'relative', width: `${width}px`, height: `${height}px`,
    overflow: 'hidden', cursor: 'crosshair'}});
  const dialog = ui.dialog({title: 'Molecule'}).add(view).showModal(true);

  const dpr = window.devicePixelRatio || 1;
  const k = Math.min(DIALOG_RESOLUTION, Math.sqrt(DIALOG_MAX_PIXELS / (width * height * dpr * dpr)));
  // uncached: Chem's renderer would otherwise keep a canvas and a picture this large for the session (Chem's own call,
  // not grok.chem.drawMolecule: a platform's js-api older than the `cache` argument would drop it)
  const molBox = await grok.functions.call('Chem:drawMolecule',
    {molStr: mol, w: Math.round(width * k), h: Math.round(height * k), popupMenu: false, cache: false}) as HTMLElement;
  const source = molBox.querySelector('canvas');
  if (!source)
    throw new Error('Chem:drawMolecule returned no canvas');
  source.style.width = `${width}px`;
  source.style.height = `${height}px`;

  const lensW = Math.round(width * LENS_FRACTION);
  const lensH = Math.round(height * LENS_FRACTION);
  const lens = ui.canvas(Math.round(lensW * dpr), Math.round(lensH * dpr));
  Object.assign(lens.style, {position: 'absolute', width: `${lensW}px`, height: `${lensH}px`, pointerEvents: 'none',
    border: '1px solid var(--grey-4)', boxShadow: '0 2px 8px rgba(0, 0, 0, 0.25)', display: 'none'});
  view.append(molBox, lens);

  // the lens follows the pointer: the source canvas' pixels around it, scaled up (the drawing may still be under way;
  // the lens reads whatever is drawn)
  const drawLens = (e: MouseEvent) => {
    const viewRect = view.getBoundingClientRect();
    const srcRect = source.getBoundingClientRect();
    lens.style.left = `${e.clientX - viewRect.left - lensW / 2}px`;
    lens.style.top = `${e.clientY - viewRect.top - lensH / 2}px`;
    lens.style.display = '';
    const s = source.width / srcRect.width; // source pixels per shown pixel
    const sw = lensW / LENS_MAGNIFICATION * s;
    const sh = lensH / LENS_MAGNIFICATION * s;
    const g = lens.getContext('2d')!;
    g.fillStyle = 'white';
    g.fillRect(0, 0, lens.width, lens.height);
    g.imageSmoothingQuality = 'high';
    g.drawImage(source, (e.clientX - srcRect.left) * s - sw / 2, (e.clientY - srcRect.top) * s - sh / 2, sw, sh,
      0, 0, lens.width, lens.height);
  };
  let frame = 0;
  view.addEventListener('mousemove', (e) => {
    cancelAnimationFrame(frame);
    frame = requestAnimationFrame(() => drawLens(e));
  });
  view.addEventListener('mouseleave', () => {
    cancelAnimationFrame(frame);
    lens.style.display = 'none';
  });
  return dialog;
}

const SCROLL_WINDOW_SIZE = 30;
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
      // Chem's one entry point for 3D coordinates (it follows Chem's engine setting)
      molBlock3D = (await grok.functions.call('Chem:embed3d', {molecule: molBlock2D})) as unknown as string;
      // NGL reads V2000 (sdf); the assembler's module converts it (no OCL)
      try {
        const assembler = await getMacromoleculeAssembler();
        molBlock3D = await assembler.toV2000(molBlock3D);
      } catch (err: any) {
        console.warn('Failed to convert molv3 to molv2', err);
      }

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
