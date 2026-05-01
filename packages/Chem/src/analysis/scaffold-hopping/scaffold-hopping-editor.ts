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

/** Parameters captured by the Scaffold Hopping dialog. */
export type ScaffoldHoppingParams = {
  table: DG.DataFrame;
  molecules: DG.Column;
  referenceRowIdx: number;
  tanimotoMin: number;
  tanimotoMax: number;
  mcsRatioMax: number;
  minPharmOverlap: number;
  /** RDKit atom indices the user marked as the *replaceable* region — i.e.
   *  the part the scaffold hop should change. Empty means "let the system
   *  decide" (auto-Murcko fallback). */
  replaceableAtoms: number[];
  /** When true, the SH flag additionally requires the candidate's ECFP4 Tc
   *  to be inside the pre-filter window. When false, only the Maeda MCS
   *  atom-ratio (and Pharm Sim if enabled) drive the flag — closer to the
   *  paper's pure SH definition. The Tc range is still used as a performance
   *  shortlist regardless of this flag. */
  useTcInFlag: boolean;
  /** When true, the SH flag additionally requires Pharm Sim ≥ the threshold.
   *  When false, only the Maeda MCS atom-ratio (and Tc window if enabled)
   *  drive the flag. */
  usePharmInFlag: boolean;
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
 *    (atom-ratio + Pharm overlap min) sections.
 *
 *  The composite ranking score combines three descriptors automatically
 *  (`0.4 × Tc + 0.3 × Pharm2D + 0.3 × CATS2D`); only Tc and Pharm have
 *  user-tunable thresholds because they participate in the flag. CATS
 *  contributes to the ranking score only — surfaced as a `Scaffold Hop
 *  CATS Sim` column on the result table. */
export class ScaffoldHoppingFunctionEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot: HTMLElement = ui.div();
  refRowInput!: DG.InputBase<number | null>;

  tanimotoMinInput = ui.input.float('ECFP4 Tc min', {value: 0.2,
    tooltipText: 'Pre-filter shortlist lower bound. Drops rows with negligible ' +
      'similarity before pharmacophore / MCS scoring. Calibrated on Maeda 2024 ' +
      'kinase / GPCR / enzyme targets — widen for nuclear receptors or peptide ' +
      'ligands where genuine hops can fall below 0.2.'});
  tanimotoMaxInput = ui.input.float('ECFP4 Tc max', {value: 0.6,
    tooltipText: 'Pre-filter shortlist upper bound. Drops near-duplicates of ' +
      'the reference (Tc above this are usually trivial analogues, not hops).'});
  mcsMaxInput = ui.input.float('MCS atom ratio max', {value: 0.4,
    tooltipText: 'Maeda 2024 scaffold-hop criterion (J. Chem. Inf. Model. 2024, ' +
      '64, 5557): a candidate is flagged as a hop iff the atom-ratio ' +
      'ratio_atom = atoms(MCS) / atoms(reference) is ≤ this threshold. ' +
      'Default 0.4 follows the paper exactly. (Note: this is Maeda\'s SH ' +
      'classifier — distinct from TcMCS, the bond-Tanimoto formula Maeda also ' +
      'defines but uses only for chemical-space-network visualisation.)'});
  useTcInFlagInput = ui.input.bool('Tc window in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires the ' +
      'candidate\'s ECFP4 Tc to be inside [Tc min, Tc max]. When unchecked, ' +
      'the Tc range is only used as a pre-filter shortlist for performance — ' +
      'the flag itself ignores it. Uncheck this and "Pharm Sim in flag" ' +
      'together to get pure Maeda 2024 behaviour (atom-ratio ≤ 0.4 only).'});
  usePharmInFlagInput = ui.input.bool('Pharm Sim in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires Pharm Sim ' +
      '≥ Pharm Tc min. When unchecked, pharmacophore similarity is computed ' +
      'and shown but does not affect the flag — closer to Maeda\'s paper-' +
      'faithful definition (atom-ratio ≤ 0.4 only).'});

  pharmInput = ui.input.float('Pharm Tc min', {value: 0.3,
    tooltipText: 'Minimum Pharm2D (Gobbi-Poppinger 1998) pharmacophore ' +
      'Tanimoto similarity between the candidate and the reference. Built on ' +
      'the Chem package\'s 7-family SMARTS (Donor / Acceptor / Hydrophobic / ' +
      'Aromatic / Positive / Negative / Halogen Bond — same SMARTS the ' +
      'Pharmacophore Features info panel uses) with feature pairs AND triplets ' +
      '(maxPointCount=3) at distance bins (0,2)(2,4)(4,6)(6,8)(8,12). Encodes ' +
      'where features sit relative to each other on the molecular graph, not ' +
      'just whether they exist. Computed server-side via RDKit Chem.Pharm2D.'});

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
    'Remove all marked atoms — the system will fall back to auto-Murcko');

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
    this.pharmInput.addValidator(this._rangeValidator);
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
      this.selectionCaption.innerHTML = 'Click atoms to mark the replaceable region. ' +
        '<span style="color:var(--grey-4)">Ctrl-hover paints, Ctrl+Shift-hover erases.</span>';
      this.selectionCaption.style.color = 'var(--grey-5)';
      this.clearSelectionBtn.style.display = 'none';
    } else {
      this.selectionCaption.textContent = `${n} atom${n === 1 ? '' : 's'} marked as replaceable.`;
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
      this.usePharmInFlagInput,
      this.pharmInput,
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
      minPharmOverlap: this.pharmInput.value!,
      replaceableAtoms: [...this.selectedAtoms].sort((a, b) => a - b),
      useTcInFlag: this.useTcInFlagInput.value ?? true,
      usePharmInFlag: this.usePharmInFlagInput.value ?? true,
    };
  }
}
