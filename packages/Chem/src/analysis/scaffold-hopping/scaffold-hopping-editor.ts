import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getRdKitModule, RDKIT_COMMON_RENDER_OPTS} from '../../utils/chem-common-rdkit';
import {extractAtomPositionsFromSvg, findNearestAtom, computeSelectedBonds}
  from '../../utils/chem-atom-picker-utils';

/** RGB triple in [0, 1] used by RDKit's `get_svg_with_highlights` for the
 *  user-marked replaceable region. Orange — visually distinct from the default
 *  RDKit query/match colours. */
const REPLACEABLE_COLOR: [number, number, number] = [1.0, 0.6, 0.2];

/** SVG canvas size for the reference preview. Square; RDKit will fit the
 *  molecule into the box and the picker uses these as click coordinates. */
const PREVIEW_SIZE = 320;

/** Three named threshold presets keyed to the Sun-Tawa-Wallqvist 2012 hop
 *  taxonomy ("Classification of scaffold hopping approaches", DDT 17:310).
 *  Each preset sets the ECFP4 Tc range and the Maeda 2024 atom-ratio cap to
 *  match the structural-change regime users typically have in mind:
 *
 *  - **Easy** (substituent / heterocycle swap, Sun-Tawa-Wallqvist class 1°):
 *    candidates share most of the scaffold; hops are decorations or single-
 *    ring isosteres. Tc 0.4-0.7 is the typical range; mcsRatio ≤ 0.7 admits
 *    most of the molecule into the MCS.
 *  - **Middle** (ring open/close, residue swap, class 2°): meaningful but
 *    bounded scaffold change. The default, balanced for typical SAR-table
 *    workflows.
 *  - **Hard / Maeda-faithful** (large topological change, class 3°/4°):
 *    paper-faithful Maeda 2024 thresholds. The Tc lower bound is dropped to
 *    0.05 because Schneider 1999, Vogt et al. 2010, and Iktos/Pinel 2023 all
 *    show genuine large-topology hops sit in the 0.05-0.30 ECFP4 range —
 *    higher cutoffs systematically discard exactly the candidates Maeda's
 *    atom-ratio classifier is designed to find. */
type PresetKey = 'Easy' | 'Middle' | 'Hard';
const PRESETS: Record<PresetKey, {tcMin: number; tcMax: number; mcsRatioMax: number}> = {
  'Easy':   {tcMin: 0.4,  tcMax: 0.7, mcsRatioMax: 0.7},
  'Middle': {tcMin: 0.2,  tcMax: 0.5, mcsRatioMax: 0.5},
  'Hard':   {tcMin: 0.05, tcMax: 0.3, mcsRatioMax: 0.4},
};
/** Display labels for the preset dropdown — keep the taxonomy hint visible
 *  so the user doesn't have to read the tooltip to understand the choice. */
const PRESET_LABELS: Record<PresetKey, string> = {
  'Easy':   'Easy — substituent / heterocycle swap',
  'Middle': 'Middle — ring open/close, residue swap',
  'Hard':   'Hard — Maeda-faithful, large topological change',
};
const DEFAULT_PRESET: PresetKey = 'Middle';

/** Parameters captured by the Scaffold Hopping dialog. */
export type ScaffoldHoppingParams = {
  table: DG.DataFrame;
  molecules: DG.Column;
  referenceRowIdx: number;
  tanimotoMin: number;
  tanimotoMax: number;
  mcsRatioMax: number;
  minCatsSim: number;
  /** RDKit atom indices the user marked as the *replaceable* region — the
   *  part the scaffold hop should change. Acts as an additional filter on
   *  top of the Maeda atom-ratio classifier: when non-empty, a candidate is
   *  flagged as a hop only if its MCS does NOT cover all marked atoms (i.e.
   *  the candidate actually changes the marked region). Empty list = no
   *  region constraint, the pipeline returns all global scaffold hops vs.
   *  the reference — paper-faithful Maeda 2024 behaviour. */
  replaceableAtoms: number[];
  /** When true, the SH flag additionally requires the candidate's ECFP4 Tc
   *  to be inside the pre-filter window. When false, only the Maeda MCS
   *  atom-ratio (and CATS Sim if enabled) drive the flag — closer to the
   *  paper's pure SH definition. The Tc range is still used as a performance
   *  shortlist regardless of this flag. */
  useTcInFlag: boolean;
  /** When true, the SH flag additionally requires CATS Sim ≥ the threshold.
   *  When false, only the Maeda MCS atom-ratio (and Tc window if enabled)
   *  drive the flag. */
  useCatsInFlag: boolean;
};

/** Editor for `Chem | Analyze | Scaffold Hopping...`.
 *
 *  Builds the popup UI:
 *  - Reference molecule preview rendered as inline RDKit SVG with paint-on-
 *    hover atom selection. The user marks the *replaceable* region (the part
 *    they want the hop to change) by hovering with Ctrl held; bonds whose
 *    endpoints are both selected are coloured automatically.
 *  - Screening-library column picker (semType: Molecule) and reference-row
 *    index input — defaults to the dataframe's `currentRowIdx`.
 *  - Pre-filter shortlist (ECFP4 Tc min/max) and Maeda 2024 hop criterion
 *    (atom-ratio + CATS Sim min) sections.
 *
 *  The composite ranking score combines two complementary descriptors
 *  automatically (`0.4 × Tc + 0.6 × CATS2D`). ECFP4 captures scaffold/
 *  topology, CATS captures the topological-pharmacophore-pair distribution
 *  that survives scaffold swaps (Schneider 1999) — the canonical low-Tc
 *  scaffold-hop signal. Both Tc and CATS Sim have user-tunable thresholds
 *  that drive the flag. */
export class ScaffoldHoppingFunctionEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot: HTMLElement = ui.div();
  refRowInput!: DG.InputBase<number | null>;

  presetInput = ui.input.choice<string>('Preset', {
    value: PRESET_LABELS[DEFAULT_PRESET],
    items: (Object.keys(PRESETS) as PresetKey[]).map((k) => PRESET_LABELS[k]),
    onValueChanged: () => this._applyPreset(),
    tooltipText: 'Picks Tc range and MCS atom-ratio cap to match the kind of ' +
      'scaffold change you\'re looking for, anchored to the Sun-Tawa-Wallqvist ' +
      '2012 taxonomy (DDT 17:310). Pick Easy for substituent / heterocycle ' +
      'swaps, Middle for ring open/close or residue swaps, Hard for large ' +
      'topological changes (paper-faithful Maeda 2024). The three numeric ' +
      'inputs below stay editable — the preset just seeds them.',
  });

  tanimotoMinInput = ui.input.float('ECFP4 Tc min', {value: PRESETS[DEFAULT_PRESET].tcMin,
    tooltipText: 'Pre-filter shortlist lower bound. Drops rows with negligible ' +
      'similarity before pharmacophore / MCS scoring. Lower this if your ' +
      'reference and the genuine hops are expected to share little topology — ' +
      'Schneider 1999 and Iktos/Pinel 2023 both report real hops at Tc < 0.2.'});
  tanimotoMaxInput = ui.input.float('ECFP4 Tc max', {value: PRESETS[DEFAULT_PRESET].tcMax,
    tooltipText: 'Pre-filter shortlist upper bound. Drops near-duplicates of ' +
      'the reference (Tc above this are usually trivial analogues, not hops).'});
  mcsMaxInput = ui.input.float('MCS atom ratio max', {value: PRESETS[DEFAULT_PRESET].mcsRatioMax,
    tooltipText: 'Maeda 2024 scaffold-hop criterion (J. Chem. Inf. Model. 2024, ' +
      '64, 5557): a candidate is flagged as a hop iff the atom-ratio ' +
      'ratio_atom = atoms(MCS) / atoms(reference) is ≤ this threshold. ' +
      'Hard preset sets 0.4 (paper-faithful); looser presets relax it.'});
  useTcInFlagInput = ui.input.bool('Tc window in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires the ' +
      'candidate\'s ECFP4 Tc to be inside [Tc min, Tc max]. When unchecked, ' +
      'the Tc range is only used as a pre-filter shortlist for performance — ' +
      'the flag itself ignores it. Uncheck this and "CATS Sim in flag" ' +
      'together to get pure Maeda 2024 behaviour (atom-ratio ≤ 0.4 only).'});
  useCatsInFlagInput = ui.input.bool('CATS Sim in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires CATS Sim ' +
      '≥ CATS Sim min. When unchecked, CATS similarity is computed and shown ' +
      'but does not affect the flag — closer to Maeda\'s paper-faithful ' +
      'definition (atom-ratio ≤ 0.4 only).'});

  catsSimInput = ui.input.float('CATS Sim min', {value: 0.8,
    tooltipText: 'Minimum CATS2D (Schneider 1999) topological-pharmacophore-' +
      'pair cosine similarity between the candidate and the reference. CATS2D ' +
      'is a 7×7×10 = 490-dim float vector counting (familyA, familyB, distance) ' +
      'pair occurrences over the Chem package\'s 7-family SMARTS (Donor / ' +
      'Acceptor / Hydrophobic / Aromatic / Positive / Negative / Halogen Bond), ' +
      'normalised by Schneider scaling so it stays scale-invariant across ' +
      'molecule sizes. Designed specifically for low-Tc scaffold-hop retrieval — ' +
      'two molecules with the same pharmacophore arrangement on different ' +
      'scaffolds will score high here even when ECFP4 says they are far apart. ' +
      'Cosine values typically run 0.7+ for related chemotypes; 0.8 default ' +
      'is permissive enough to catch genuine hops without admitting noise.'});

  /** Outer container holding the rendered SVG + click overlay. */
  referencePreview: HTMLElement = ui.div([], {style: {
    width: `${PREVIEW_SIZE}px`, height: `${PREVIEW_SIZE}px`,
    display: 'flex', alignItems: 'center', justifyContent: 'center',
    border: '1px solid var(--grey-2)', borderRadius: '4px',
    background: 'var(--white)', cursor: 'pointer', overflow: 'hidden',
  }});
  referenceCaption: HTMLElement = ui.divText('', {style: {
    fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px',
    maxWidth: `${PREVIEW_SIZE}px`, wordBreak: 'break-all',
  }});
  selectionCaption: HTMLElement = ui.divText('', {style: {
    fontSize: '11px', color: 'var(--grey-6)', marginTop: '4px', fontWeight: '500',
    maxWidth: `${PREVIEW_SIZE}px`,
  }});
  clearSelectionBtn = ui.button('Clear selection', () => this._clearSelection(),
    'Remove all marked atoms — the run will return any global scaffold hop');

  /** Atom indices the user has marked. Reset whenever the reference changes. */
  selectedAtoms: Set<number> = new Set();

  private _atomPositions: Map<number, {x: number; y: number}> = new Map();
  private _bondAtoms: Map<number, [number, number]> = new Map();
  private _svgEl: SVGSVGElement | null = null;
  private _currentRefSmiles: string | null = null;

  /** Last atom acted on by hover-with-modifier — dedup guard so the same atom
   *  isn't re-painted on every pixel of mousemove. Mirrors `_lastHoveredAtom`
   *  in `atom-picker-controller.ts`. */
  private _lastPaintedAtom: number | null = null;
  private _lastPaintMode: 'add' | 'erase' | null = null;

  constructor() {
    this.tableInput = ui.input.table('Table', {
      value: grok.shell.tv?.dataFrame ?? null,
      items: grok.shell.tables,
      onValueChanged: () => this.onTableInputChanged(),
    });
    this.tanimotoMinInput.addValidator(this._rangeValidator);
    this.tanimotoMaxInput.addValidator(this._rangeValidator);
    this.mcsMaxInput.addValidator(this._rangeValidator);
    this.catsSimInput.addValidator(this._rangeValidator);
    this.clearSelectionBtn.style.display = 'none';
    this.onTableInputChanged();
  }

  private _rangeValidator = (s: string): string | null => {
    const v = parseFloat(s);
    if (Number.isNaN(v)) return 'Number expected';
    if (v < 0 || v > 1) return 'Must be in [0, 1]';
    return null;
  };

  onTableInputChanged() {
    const table = this.tableInput.value;
    if (!table) return;

    const molColumns = table.columns.toList().filter((c) => c.semType === DG.SEMTYPE.MOLECULE);
    const firstMol = molColumns[0];

    const newColInput = ui.input.column('Screening library', {
      table,
      value: firstMol,
      filter: (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE,
      onValueChanged: () => this._refreshReferencePreview(),
    });
    ui.empty(this.colInputRoot);
    Array.from(newColInput.root.children).forEach((c) => this.colInputRoot.append(c));
    this.colInput = newColInput;

    const defaultIdx = table.currentRowIdx >= 0 ? table.currentRowIdx : 0;
    this.refRowInput = ui.input.int('Reference row', {
      value: defaultIdx,
      onValueChanged: () => this._refreshReferencePreview(),
    });
    this.refRowInput.addValidator((s: string) => {
      const v = parseInt(s, 10);
      if (Number.isNaN(v)) return 'Integer expected';
      if (v < 0 || v >= table.rowCount) return `Must be in [0, ${table.rowCount - 1}]`;
      return null;
    });

    this._refreshReferencePreview();
  }

  private _refreshReferencePreview() {
    const table = this.tableInput.value;
    const col = this.colInput?.value;
    const rowIdx = this.refRowInput?.value;
    this.selectedAtoms.clear();
    this._lastPaintedAtom = null;
    this._lastPaintMode = null;
    this._atomPositions.clear();
    this._bondAtoms.clear();
    this._svgEl = null;
    ui.empty(this.referencePreview);

    if (!table || !col || rowIdx == null || rowIdx < 0 || rowIdx >= table.rowCount) {
      this.referencePreview.appendChild(ui.divText('No reference molecule',
        {style: {color: 'var(--grey-4)'}}));
      this.referenceCaption.textContent = '';
      this._currentRefSmiles = null;
      this._updateSelectionCaption();
      return;
    }
    const smiles = col.get(rowIdx);
    this._currentRefSmiles = smiles;
    if (!smiles) {
      this.referencePreview.appendChild(ui.divText('Empty cell',
        {style: {color: 'var(--grey-4)'}}));
      this.referenceCaption.textContent = `Row ${rowIdx} • (empty)`;
      this._updateSelectionCaption();
      return;
    }
    this.referenceCaption.textContent = `Row ${rowIdx} • ${smiles}`;
    this._renderInteractivePicker(smiles);
  }

  private _renderInteractivePicker(smiles: string) {
    let mol: any = null;
    try {
      const rdKit = getRdKitModule();
      mol = rdKit.get_mol(smiles);
      const svgString = this._buildSvgWithHighlights(mol);
      const doc = new DOMParser().parseFromString(svgString, 'image/svg+xml');
      const svgEl = doc.documentElement as unknown as SVGSVGElement;
      this._svgEl = svgEl;
      svgEl.setAttribute('width', `${PREVIEW_SIZE}`);
      svgEl.setAttribute('height', `${PREVIEW_SIZE}`);
      svgEl.style.cursor = 'pointer';

      this.referencePreview.appendChild(svgEl);

      const {positions, bondAtoms} = extractAtomPositionsFromSvg(svgEl);
      this._atomPositions = positions;
      this._bondAtoms = bondAtoms;

      this._attachPickerListeners(svgEl, smiles);
    } catch (e: any) {
      this.referencePreview.appendChild(ui.divText(`Render error: ${e?.message ?? e}`,
        {style: {color: 'var(--failure)', fontSize: '11px', padding: '8px'}}));
    } finally {
      mol?.delete();
    }
    this._updateSelectionCaption();
  }

  /** Renders the molecule with `RDKIT_COMMON_RENDER_OPTS` so it matches the
   *  grid cell renderer's appearance (line widths, font sizes, atom palette). */
  private _buildSvgWithHighlights(mol: any): string {
    const atomIdxs = [...this.selectedAtoms];
    const highlightAtomColors: {[k: number]: number[]} = {};
    for (const a of atomIdxs) highlightAtomColors[a] = [...REPLACEABLE_COLOR];
    const {bondsArr, highlightBondColors} = computeSelectedBonds(
      this.selectedAtoms, this._bondAtoms, [...REPLACEABLE_COLOR]);
    const details: {[key: string]: any} = {
      ...RDKIT_COMMON_RENDER_OPTS,
      width: PREVIEW_SIZE,
      height: PREVIEW_SIZE,
      atoms: atomIdxs,
      bonds: bondsArr,
      highlightAtomColors,
      highlightBondColors,
    };
    return mol.get_svg_with_highlights(JSON.stringify(details));
  }

  private _attachPickerListeners(svgEl: SVGSVGElement, smiles: string) {
    svgEl.addEventListener('click', (ev: MouseEvent) => this._onSvgClick(ev, smiles));
    svgEl.addEventListener('mousemove', (ev: MouseEvent) => this._onSvgMouseMove(ev, smiles));
    svgEl.addEventListener('mouseleave', () => this._resetPaintTracking());
  }

  private _svgClientToLocal(ev: MouseEvent): {x: number; y: number} | null {
    if (!this._svgEl) return null;
    const ctm = this._svgEl.getScreenCTM();
    if (!ctm) return null;
    const pt = this._svgEl.createSVGPoint();
    pt.x = ev.clientX;
    pt.y = ev.clientY;
    const local = pt.matrixTransform(ctm.inverse());
    return {x: local.x, y: local.y};
  }

  /** Click handler — modifier semantics mirror `atom-picker-controller.ts:203-204`:
   *  - **Ctrl/Cmd + click**: add the atom (idempotent — no-op if already in)
   *  - **Ctrl/Cmd + Shift + click**: remove the atom (idempotent)
   *  - **Plain click**: toggle */
  private _onSvgClick(ev: MouseEvent, smiles: string) {
    if (!this._svgEl || smiles !== this._currentRefSmiles) return;
    ev.preventDefault();
    const local = this._svgClientToLocal(ev);
    if (!local) return;
    const nearest = findNearestAtom(this._atomPositions, local.x, local.y);
    if (nearest === null) return;

    const isErase = (ev.ctrlKey || ev.metaKey) && ev.shiftKey;
    const isPaint = (ev.ctrlKey || ev.metaKey) && !ev.shiftKey;

    let changed = false;
    if (isErase) {
      changed = this.selectedAtoms.delete(nearest);
    } else if (isPaint) {
      if (!this.selectedAtoms.has(nearest)) {
        this.selectedAtoms.add(nearest);
        changed = true;
      }
    } else {
      if (this.selectedAtoms.has(nearest)) this.selectedAtoms.delete(nearest);
      else this.selectedAtoms.add(nearest);
      changed = true;
    }
    if (changed) this._reRender(smiles);
  }

  /** Hover handler — paints atoms directly into the selection while the
   *  modifier is held, the same as `rdkit-cell-renderer.ts:723-726`:
   *  Ctrl/Cmd+hover adds, Ctrl/Cmd+Shift+hover removes, no modifier no-ops.
   *  Dedup via `_lastPaintedAtom` so the cursor wiggling inside one atom's
   *  hit zone doesn't fire repeated re-renders. */
  private _onSvgMouseMove(ev: MouseEvent, smiles: string) {
    if (smiles !== this._currentRefSmiles) return;
    const isErase = (ev.ctrlKey || ev.metaKey) && ev.shiftKey;
    const isPaint = (ev.ctrlKey || ev.metaKey) && !ev.shiftKey;

    if (!isErase && !isPaint) {
      this._lastPaintedAtom = null;
      this._lastPaintMode = null;
      return;
    }
    const local = this._svgClientToLocal(ev);
    if (!local) return;
    const nearest = findNearestAtom(this._atomPositions, local.x, local.y);
    if (nearest === null) {
      this._lastPaintedAtom = null;
      return;
    }
    const mode: 'add' | 'erase' = isErase ? 'erase' : 'add';
    if (nearest === this._lastPaintedAtom && mode === this._lastPaintMode) return;
    this._lastPaintedAtom = nearest;
    this._lastPaintMode = mode;

    let changed = false;
    if (mode === 'erase') {
      changed = this.selectedAtoms.delete(nearest);
    } else {
      if (!this.selectedAtoms.has(nearest)) {
        this.selectedAtoms.add(nearest);
        changed = true;
      }
    }
    if (changed) this._reRender(smiles);
  }

  private _resetPaintTracking() {
    this._lastPaintedAtom = null;
    this._lastPaintMode = null;
  }

  /** Applies the currently-selected preset's Tc/mcsRatio values to the three
   *  numeric inputs. The user can still tweak the inputs manually after the
   *  preset is applied — we don't lock them. */
  private _applyPreset() {
    const label = this.presetInput.value;
    const key = (Object.keys(PRESETS) as PresetKey[]).find((k) => PRESET_LABELS[k] === label);
    if (!key) return;
    const p = PRESETS[key];
    this.tanimotoMinInput.value = p.tcMin;
    this.tanimotoMaxInput.value = p.tcMax;
    this.mcsMaxInput.value = p.mcsRatioMax;
  }

  private _reRender(smiles: string) {
    let mol: any = null;
    try {
      mol = getRdKitModule().get_mol(smiles);
      const svgString = this._buildSvgWithHighlights(mol);
      const doc = new DOMParser().parseFromString(svgString, 'image/svg+xml');
      const newSvg = doc.documentElement as unknown as SVGSVGElement;
      newSvg.setAttribute('width', `${PREVIEW_SIZE}`);
      newSvg.setAttribute('height', `${PREVIEW_SIZE}`);
      newSvg.style.cursor = 'pointer';
      this._attachPickerListeners(newSvg, smiles);
      ui.empty(this.referencePreview);
      this.referencePreview.appendChild(newSvg);
      this._svgEl = newSvg;
    } catch (_) {
      // leave the previous render in place
    } finally {
      mol?.delete();
    }
    this._updateSelectionCaption();
  }

  private _clearSelection() {
    if (this.selectedAtoms.size === 0) return;
    this.selectedAtoms.clear();
    if (this._currentRefSmiles) this._reRender(this._currentRefSmiles);
    else this._updateSelectionCaption();
  }

  private _updateSelectionCaption() {
    const n = this.selectedAtoms.size;
    if (n === 0) {
      this.selectionCaption.innerHTML =
        '<span style="color:var(--grey-6)">Optional: mark a region to find hops that change it.</span><br>' +
        '<span style="color:var(--grey-4)">Ctrl-hover paints, Ctrl+Shift-hover erases. ' +
        'Leave empty to find any global scaffold hop.</span>';
      this.selectionCaption.style.color = '';
      this.clearSelectionBtn.style.display = 'none';
    } else {
      this.selectionCaption.textContent =
        `${n} atom${n === 1 ? '' : 's'} marked — hits will be filtered to hops that change this region.`;
      this.selectionCaption.style.color = 'var(--orange-2, #c87325)';
      this.clearSelectionBtn.style.display = '';
    }
  }

  public getEditor(): HTMLElement {
    const left = ui.divV([
      this.referencePreview,
      this.referenceCaption,
      this.selectionCaption,
      this.clearSelectionBtn,
    ], {style: {paddingRight: '12px'}});

    const right = ui.divV([
      this.tableInput,
      this.colInputRoot,
      this.refRowInput,
      ui.h3('Threshold preset',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.presetInput,
      ui.divText('Sets Tc range + MCS atom-ratio cap. The numeric inputs ' +
        'below remain editable; the preset just seeds them. Anchored to the ' +
        'Sun-Tawa-Wallqvist 2012 hop taxonomy (DDT 17:310).',
        {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),
      ui.h3('Pre-filter shortlist (ECFP4 Tanimoto)',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.tanimotoMinInput,
      this.tanimotoMaxInput,
      ui.h3('Hop criterion (Maeda 2024 — always applied)',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.mcsMaxInput,
      ui.divText('Maeda 2024 SH classifier (J. Chem. Inf. Model. 64, 5557): ' +
        'ratio_atom = atoms(MCS) / atoms(reference). A row passes this ' +
        'criterion iff ratio_atom ≤ max.',
        {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),

      ui.h3('Optional flag conditions',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.useCatsInFlagInput,
      this.catsSimInput,
      this.useTcInFlagInput,
      ui.divText('Each unchecked condition is computed and shown but does not ' +
        'affect the Scaffold Hop boolean flag. Uncheck both to reproduce the ' +
        'Maeda 2024 paper\'s definition (atom-ratio ≤ 0.4 only).',
        {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),
    ], {style: {minWidth: '340px'}, classes: 'ui-form'});

    return ui.divH([left, right], {style: {alignItems: 'flex-start'}});
  }

  public getParams(): ScaffoldHoppingParams {
    const table = this.tableInput.value!;
    const molecules = this.colInput.value!;
    const referenceRowIdx = this.refRowInput.value!;
    if (this.tanimotoMinInput.value! > this.tanimotoMaxInput.value!)
      throw new Error('ECFP4 Tc min must be ≤ ECFP4 Tc max');
    return {
      table, molecules, referenceRowIdx,
      tanimotoMin: this.tanimotoMinInput.value!,
      tanimotoMax: this.tanimotoMaxInput.value!,
      mcsRatioMax: this.mcsMaxInput.value!,
      minCatsSim: this.catsSimInput.value!,
      replaceableAtoms: [...this.selectedAtoms].sort((a, b) => a - b),
      useTcInFlag: this.useTcInFlagInput.value ?? true,
      useCatsInFlag: this.useCatsInFlagInput.value ?? true,
    };
  }
}
