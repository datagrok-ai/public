/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import $ from 'cash-dom';
import * as rxjs from 'rxjs';

import * as C from '../utils/constants';
import {PeptideUtils} from '../peptideUtils';
import {scaleActivity, highlightMonomerPosition} from '../utils/misc';
import {calculateMonomerPositionStatistics} from '../utils/algorithms';
import {MonomerPositionStats, PositionStats, StatsItem} from '../utils/statistics';
import {showTooltip} from '../utils/tooltips';
import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {StringListSeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {ALPHABET, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GapOriginals, GAP_SYMBOL, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {SARViewer} from './sar-viewer';

export enum PEPTIDE_GENERATION_PROPS {
  SEQUENCE = 'sequence',
  ACTIVITY = 'activity',
  ACTIVITY_SCALING = 'activityScaling',
  ACTIVITY_TARGET = 'activityTarget',
  PEPTIDE_COUNT = 'peptideCount',
  CANDIDATES_PER_POSITION = 'candidatesPerPosition',
  BEAM_WIDTH = 'beamWidth',
  MIN_SUPPORT = 'minSupport',
  MAX_P_VALUE = 'maxPValue',
  SHRINKAGE = 'shrinkageStrength',
  EXCLUDE_EXISTING = 'excludeExisting',
  ONLY_SIGNIFICANT = 'onlySignificantPositions',
}

/** Generated-peptides result grid column names. */
export enum GEN_COLUMN_NAMES {
  PEPTIDE = 'Peptide',
  PREDICTED = 'Predicted activity',
  CONFIDENCE = 'Confidence',
  MIN_SUPPORT = 'Min support',
  MEAN_P_VALUE = 'Mean p-value',
  SIGNIFICANT_POSITIONS = 'Significant positions',
  NOVEL = 'Novel',
  BASIS = 'Basis',
}

/** A single candidate monomer choice at a given position. */
type Candidate = {
  position: string,
  monomer: string,
  /** Deviation of this monomer's group mean activity from the global mean (additive contribution). */
  contribution: number,
  count: number,
  pValue: number | null,
  significant: boolean,
  supported: boolean,
};

/** A partial peptide during beam search. */
type BeamItem = {
  monomers: string[],
  candidates: Candidate[],
  cumContribution: number,
};

export const PEPTIDE_GENERATION_VIEWER_NAME = 'Peptide Generation';

/**
 * De novo peptide generator.
 *
 * Reuses the per-monomer-per-position statistics (the same engine that powers the Most Potent Residues and
 * Sequence Variability Map viewers) to synthesize a ranked set of new candidate peptides under an additive
 * positional (Free-Wilson-style) activity model:
 *
 *   predicted(peptide) = globalMean + Σ_positions (meanActivity(monomer@position) − globalMean)
 *
 * Candidates are enumerated with a bounded beam search over the best-performing monomers at each position and
 * ranked by predicted activity for the chosen target (High/Low). Each generated peptide is reported together
 * with the statistical support behind the prediction, so the (position-independence) assumption of the model
 * stays visible to the user.
 */
export class PeptideGenerationViewer extends DG.JsViewer {
  sequenceColumnName: string;
  activityColumnName: string;
  activityScaling: C.SCALING_METHODS;
  activityTarget: C.ACTIVITY_TARGET;
  peptideCount: number;
  candidatesPerPosition: number;
  beamWidth: number;
  minSupport: number;
  maxPValue: number;
  shrinkageStrength: number;
  excludeExisting: boolean;
  onlySignificantPositions: boolean;

  private _grid: DG.Grid | null = null;
  private _positionColumns: DG.Column<string>[] | null = null;
  private _monomerPositionStats: MonomerPositionStats | null = null;
  /** Scaled activity column aligned to the source dataframe, used for the invariant-map style tooltips. */
  private _scaledActivityCol: DG.Column<number> | null = null;
  /** Per generated position (grid column name) → per result-row chosen monomer stats, for cell rendering/tooltips. */
  private _basisByPosition: {[position: string]: (Candidate | null)[]} = {};
  private _basisPositionNames: string[] = [];
  private _basisPositionSet: Set<string> = new Set();
  /** Per result-row max absolute contribution, for fast per-row circle-size normalization. */
  private _rowMaxAbsContribution: number[] = [];
  /** Monomer symbol → text color derived from the monomer library background color. */
  private _monomerColorCache: Map<string, string> = new Map();
  private _biotype: HelmType = HelmTypes.AA;

  constructor() {
    super();
    this.sequenceColumnName = this.column(PEPTIDE_GENERATION_PROPS.SEQUENCE,
      {semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.activityColumnName = this.column(PEPTIDE_GENERATION_PROPS.ACTIVITY,
      {columnTypeFilter: 'numerical', nullable: false});
    this.activityScaling = this.string(PEPTIDE_GENERATION_PROPS.ACTIVITY_SCALING, C.SCALING_METHODS.NONE,
      {choices: Object.values(C.SCALING_METHODS), nullable: false,
        description: 'Activity transformation applied before computing statistics'}) as C.SCALING_METHODS;
    this.activityTarget = this.string(PEPTIDE_GENERATION_PROPS.ACTIVITY_TARGET, C.ACTIVITY_TARGET.HIGH,
      {choices: Object.values(C.ACTIVITY_TARGET), nullable: false,
        description: 'Whether to design peptides that maximize (High) or minimize (Low) the scaled activity'}) as C.ACTIVITY_TARGET;
    this.peptideCount = this.int(PEPTIDE_GENERATION_PROPS.PEPTIDE_COUNT, 50,
      {min: 1, max: 1000, description: 'How many peptides to generate'});
    this.candidatesPerPosition = this.int(PEPTIDE_GENERATION_PROPS.CANDIDATES_PER_POSITION, 3,
      {min: 1, max: 10, description: 'Number of best-performing monomers considered at each position'});
    this.beamWidth = this.int(PEPTIDE_GENERATION_PROPS.BEAM_WIDTH, 100,
      {min: 1, max: 2000, description: 'Beam search width. Larger values explore more combinations at higher cost'});
    this.minSupport = this.int(PEPTIDE_GENERATION_PROPS.MIN_SUPPORT, 3,
      {min: 1, max: 1000, description: 'Minimum number of observed sequences backing a monomer for it to be trusted'});
    this.maxPValue = this.float(PEPTIDE_GENERATION_PROPS.MAX_P_VALUE, 1,
      {min: 0, max: 1, description: 'Only consider monomers with a t-test p-value at or below this value (1 = no filtering)'});
    this.shrinkageStrength = this.float(PEPTIDE_GENERATION_PROPS.SHRINKAGE, 5,
      {min: 0, max: 100, description: 'Shrinks each position\'s contribution towards zero by a factor of ' +
        'count/(count+k), damping poorly-supported monomers so they cannot inflate the predicted activity. ' +
        '0 disables shrinkage (raw additive model, prone to over-prediction)'});
    this.excludeExisting = this.bool(PEPTIDE_GENERATION_PROPS.EXCLUDE_EXISTING, true,
      {description: 'Exclude generated peptides that already exist in the source data'});
    this.onlySignificantPositions = this.bool(PEPTIDE_GENERATION_PROPS.ONLY_SIGNIFICANT, false,
      {description: 'Vary only positions that have at least one statistically significant monomer; keep the ' +
        'most common monomer elsewhere'});
  }

  get name(): string {
    return PEPTIDE_GENERATION_VIEWER_NAME;
  }

  onTableAttached(): void {
    super.onTableAttached();
    const seqCol = this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
    if (seqCol == null) {
      grok.shell.warning('Peptide Generation: no Macromolecule column found');
      return;
    }
    this.getProperty(PEPTIDE_GENERATION_PROPS.SEQUENCE + 'ColumnName')?.set(this, seqCol.name);
    const potentialActivity = wu(this.dataFrame.columns.numerical)
      .find((col) => col.name.toLowerCase().includes('activity')) ?? wu(this.dataFrame.columns.numerical).next()?.value;
    if (potentialActivity)
      this.getProperty(PEPTIDE_GENERATION_PROPS.ACTIVITY + 'ColumnName')?.set(this, potentialActivity.name);
    this.subs.push(DG.debounce(this.dataFrame.onFilterChanged, 500).subscribe(() => {
      this._propChangeTimer && clearTimeout(this._propChangeTimer);
      this._monomerPositionStats = null;
      this._grid = null;
      this.render();
    }));
    this.render();
  }

  private _propChangeTimer: ReturnType<typeof setTimeout> | null = null;
  onPropertyChanged(_property: DG.Property | null): void {
    if (!this.dataFrame || !this.sequenceColumnName || !this.activityColumnName)
      return;
    this._propChangeTimer && clearTimeout(this._propChangeTimer);
    this._propChangeTimer = setTimeout(() => {
      this._positionColumns = null;
      this._monomerPositionStats = null;
      this._grid = null;
      this.render();
    }, 500);
  }

  /** PeptidesModel already attached to the dataframe, or null. Never force-creates one. */
  private get existingModel(): PeptidesModel | null {
    return (this.dataFrame?.temp?.[PeptidesModel.modelName] as PeptidesModel) ?? null;
  }

  /**
   * Per-monomer-per-position statistics. Reuses the cached statistics of the SAR viewers (or the analysis model)
   * when they cover the same sequence/activity/scaling and filtering state, mirroring how those viewers share stats
   * among themselves. Falls back to computing its own (filter-aware) statistics when nothing reusable is present.
   */
  get monomerPositionStats(): MonomerPositionStats {
    if (this._monomerPositionStats != null)
      return this._monomerPositionStats;

    const wantFiltered = this.dataFrame.filter.anyFalse;
    const paramsMatch = (o: {sequenceColumnName?: string, activityColumnName?: string,
      activityScaling?: C.SCALING_METHODS} | null | undefined): boolean =>
      o != null && this.sequenceColumnName === o.sequenceColumnName &&
      this.activityColumnName === o.activityColumnName && this.activityScaling === o.activityScaling;

    const model = this.existingModel;
    if (model != null) {
      // A sibling SAR viewer may already hold matching stats. Its mean/count/p-value are pure activity statistics
      // (unaffected by any invariant-map value aggregation), so only the data and filtering intent must match.
      for (const viewerType of [VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP, VIEWER_TYPE.MOST_POTENT_RESIDUES]) {
        const viewer = model.findViewer(viewerType) as SARViewer | null;
        if (viewer?._monomerPositionStats != null && paramsMatch(viewer) &&
          (!wantFiltered || viewer.dataSource === 'Filtered'))
          return this._monomerPositionStats = viewer._monomerPositionStats;
      }
      // The model caches unfiltered statistics; safe to reuse only when we are not filtering.
      if (!wantFiltered && model.monomerPositionStats != null && paramsMatch(model.settings))
        return this._monomerPositionStats = model.monomerPositionStats;
    }

    const scaledActivityCol = scaleActivity(this.dataFrame.getCol(this.activityColumnName), this.activityScaling);
    this._monomerPositionStats = calculateMonomerPositionStatistics(scaledActivityCol, this.dataFrame.filter,
      this.getPositionColumns(), {isFiltered: wantFiltered});
    return this._monomerPositionStats;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
    this._currentCellSub?.unsubscribe();
  }

  /** Resolves per-position monomer columns, reusing the ones produced by an existing analysis when possible. */
  private getPositionColumns(): DG.Column<string>[] {
    if (this._positionColumns != null)
      return this._positionColumns;
    const tagged = wu(this.dataFrame.columns.toList())
      .filter((col) => col.getTag(C.TAGS.POSITION_COL) === `${true}`).toArray();
    if (tagged.length !== 0)
      this._positionColumns = tagged as DG.Column<string>[];
    else {
      this._positionColumns = splitAlignedSequences(this.dataFrame.getCol(this.sequenceColumnName),
        PeptideUtils.getSeqHelper()).columns.toList();
    }
    return this._positionColumns;
  }

  render(): void {
    $(this.root).empty();
    if (!this.dataFrame || !this.sequenceColumnName || !this.activityColumnName) {
      this.root.appendChild(ui.divText('Please select a sequence and an activity column in the viewer properties'));
      return;
    }

    try {
      this._grid ??= this.createGrid();
    } catch (e) {
      grok.log.error(e as string);
      this.root.appendChild(ui.divText(`Peptide Generation failed: ${(e as Error)?.message ?? e}`));
      return;
    }

    const targetLabel = this.activityTarget === C.ACTIVITY_TARGET.HIGH ? 'highest' : 'lowest';
    const subtitle = ui.divText(
      `${this._grid.table.rowCount} peptides generated for ${targetLabel} activity ` +
      `(additive positional model — treats positions as independent)`,
      {style: {color: 'var(--grey-5)', fontSize: '12px'}});

    const addIcon = ui.iconFA('plus', () => this.openAsTableView(),
      'Open the generated peptides as a new table');
    addIcon.style.cursor = 'pointer';
    addIcon.style.color = 'var(--blue-1)';

    const infoIcon = ui.iconFA('info-circle');
    infoIcon.style.color = 'var(--blue-1)';
    infoIcon.style.cursor = 'help';
    ui.tooltip.bind(infoIcon, () => this.buildInfoTooltip());

    const header = ui.divH([subtitle, addIcon, infoIcon],
      {style: {alignItems: 'center', gap: '6px', padding: '4px 8px'}});

    const gridRoot = this._grid.root;
    gridRoot.style.width = '100%';
    this.root.appendChild(ui.divV([header, gridRoot], {style: {height: '100%'}}));
    this._grid.invalidate();
  }

  /**
   * Opens the generated peptides as a standalone table view. The custom-rendered per-position basis columns are
   * flattened to plain "monomer, contribution, p-value" text so the exported table is self-describing.
   */
  private openAsTableView(): void {
    if (this._grid == null)
      return;
    const clone = this._grid.table.clone();
    clone.name = 'Generated Peptides';
    for (const posName of this._basisPositionNames) {
      const col = clone.col(posName);
      const cellData = this._basisByPosition[posName];
      if (col == null || cellData == null)
        continue;
      col.init((i) => {
        const c = cellData[i];
        if (c == null || !c.monomer)
          return '';
        const pValue = c.pValue == null ? 'n/a' :
          (c.pValue < 0.001 ? c.pValue.toExponential(1) : c.pValue.toFixed(3));
        return `${c.monomer}, ${this.formatContribution(c.contribution)}, ${pValue}`;
      });
    }
    const v = grok.shell.addTableView(clone);
    setTimeout(() => grok.shell.v = v, 50);
  }

  /** Rich explanation shown on hovering the info icon: how the algorithm works and how to read the results. */
  private buildInfoTooltip(): HTMLElement {
    return ui.divV([
      ui.h3('How Peptide Generation works'),
      ui.divText('Additive positional (Free-Wilson) model', {style: {fontWeight: 'bold'}}),
      ui.divText('Predicted activity = global mean + Σ over positions of (chosen monomer\'s mean activity − ' +
        'global mean). Each per-position contribution is shrunk by count / (count + k) so that poorly-supported ' +
        'monomers cannot dominate the prediction.'),
      ui.divText('Generation', {style: {fontWeight: 'bold', marginTop: '6px'}}),
      ui.divText('At each position the best-scoring monomers toward the chosen High/Low target are kept, then a ' +
        'beam search combines them into whole peptides, ranked by predicted activity.'),
      ui.divText('Per-position basis columns', {style: {fontWeight: 'bold', marginTop: '6px'}}),
      ui.divText('The trailing columns (one per position) show the chosen monomer, an effect circle, and the ' +
        'signed contribution. In each circle:'),
      ui.list([
        'Size = effect magnitude (|contribution|), normalized within the row: bigger = stronger effect on activity.',
        'Color = direction and significance: red raises activity, blue lowers it, pale ≈ weak significance, ' +
          'grey = not enough data to test.',
        'Hover a cell for the full activity distribution of that monomer-at-position versus the rest.',
      ], {processNode: (a) => {a.style.maxWidth = '470px'; a.style.display = 'block';},
      }),
      ui.divText('Keep in mind', {style: {fontWeight: 'bold', marginTop: '6px'}}),
      ui.list([
        'Predicted activity is a ranking score / optimistic estimate, NOT a calibrated measurement.',
        'Because contributions add across positions, it can extrapolate well beyond the observed activity ' +
          'range — a peptide combining the best residue at every position may not exist or behave additively ' +
          '(epistasis is invisible to this model).',
        'Judge reliability with the Confidence, Min support, Significant positions and Basis columns.',
        'Raise "Shrinkage strength" or enable "Only significant positions" for more conservative magnitudes. ' +
          'On data with no real signal, few or no significant peptides should appear.',
      ], {processNode: (a) => {a.style.maxWidth = '470px'; a.style.display = 'block';}}),
    ], {style: {maxWidth: '480px', fontSize: '12px'}});
  }

  private _currentCellSub: rxjs.Subscription | null = null;
  /** Builds the generated-peptides grid. */
  private createGrid(): DG.Grid {
    const resultDf = this.generatePeptides();
    const grid = resultDf.plot.grid();
    grid.props.allowEdit = false;
    grid.props.showReadOnlyNotifications = false;
    this._currentCellSub?.unsubscribe();
    this._currentCellSub = grid.onCurrentCellChanged.subscribe((g) => {
      if (g && g.dart && g.cell.dart && (g.tableRowIndex ?? -1) > -1 && g.cell.column?.name === GEN_COLUMN_NAMES.PEPTIDE)
        grok.shell.o = DG.SemanticValue.fromGridCell(g);
    });

    const predictedCol = grid.col(GEN_COLUMN_NAMES.PREDICTED);
    if (predictedCol)
      predictedCol.format = '#.000';
    const pValCol = grid.col(GEN_COLUMN_NAMES.MEAN_P_VALUE);
    if (pValCol)
      pValCol.format = '#.000';
    const confCol = grid.col(GEN_COLUMN_NAMES.CONFIDENCE);
    if (confCol)
      confCol.format = '0.##%';

    // The per-position basis columns supersede the terse text summary; keep the data but hide the column.
    const basisTextCol = grid.col(GEN_COLUMN_NAMES.BASIS);
    if (basisTextCol)
      basisTextCol.visible = false;

    // Custom-render the per-position basis cells: monomer + effect circle + contribution number.
    const alphabet = this.dataFrame.getCol(this.sequenceColumnName).getTag(bioTAGS.alphabet) ?? ALPHABET.UN;
    this._biotype = (alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA) ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
    this._monomerColorCache.clear();
    grid.props.rowHeight = 22;
    for (const posName of this._basisPositionNames) {
      const gc = grid.col(posName);
      if (gc)
        gc.width = 85;
    }
    grid.onCellRender.subscribe((args) => this.renderBasisCell(args));
    grid.onCellTooltip((gridCell, x, y) => this.showBasisTooltip(gridCell, x, y));
    grid.root.addEventListener('mouseleave', () => this.dataFrame?.rows?.highlight(() => false));

    // Highest (or lowest) predicted activity on top.
    grid.sort([GEN_COLUMN_NAMES.PREDICTED], [this.activityTarget !== C.ACTIVITY_TARGET.HIGH]);
    return grid;
  }

  /**
   * Renders a per-position basis cell: the chosen monomer, a circle whose radius encodes the (per-row normalized)
   * absolute contribution and whose color encodes direction and significance, and the signed contribution value.
   * @param args - Grid cell render arguments.
   */
  private renderBasisCell(args: DG.GridCellRenderArgs): void {
    const cell = args.cell;
    const posName = cell.tableColumn?.name;
    if (!cell.isTableCell || posName == null || !this._basisPositionSet.has(posName))
      return;

    const g = args.g;
    const b = args.bounds;
    g.save();
    g.beginPath();
    g.rect(b.x, b.y, b.width, b.height);
    g.clip();

    const rowIdx = cell.tableRowIndex;
    const cand = rowIdx == null ? null : this._basisByPosition[posName]?.[rowIdx] ?? null;
    if (cand == null || !cand.monomer) {
      args.preventDefault();
      g.restore();
      return;
    }

    const pad = 3;
    const midY = b.y + b.height / 2;
    const monomerW = b.width * 0.36;
    const numberW = b.width * 0.40;
    const circleZoneX = b.x + pad + monomerW;
    const circleZoneW = Math.max(6, b.width - monomerW - numberW - 2 * pad);

    g.textBaseline = 'middle';
    g.font = '11px Roboto, Roboto Local, sans-serif';

    // Monomer symbol, colored with its monomer-library background color, ellipsized to the allocated width.
    g.textAlign = 'left';
    g.fillStyle = this.monomerTextColor(cand.monomer);
    g.fillText(this.fitText(g, cand.monomer, monomerW - pad), b.x + pad, midY);

    // Effect circle: radius ∝ |contribution| (per-row normalized); color = direction × significance (red/blue).
    const rowMaxAbs = this._rowMaxAbsContribution[rowIdx!] || 0;
    const rCoef = rowMaxAbs > 0 ? Math.min(1, Math.abs(cand.contribution) / rowMaxAbs) : 0;
    const maxRadius = Math.min(circleZoneW, b.height) / 2 * 0.9;
    const radius = Math.max(2, maxRadius * rCoef);
    g.beginPath();
    g.fillStyle = this.contributionColor(cand.contribution, cand.pValue);
    g.arc(circleZoneX + circleZoneW / 2, midY, radius, 0, Math.PI * 2);
    g.fill();
    g.lineWidth = 0.5;
    g.strokeStyle = 'rgba(0,0,0,0.18)';
    g.stroke();

    // Signed contribution value.
    g.textAlign = 'right';
    g.fillStyle = '#606060';
    g.fillText(this.fitText(g, this.formatContribution(cand.contribution), numberW - pad), b.x + b.width - pad, midY);

    args.preventDefault();
    g.restore();
  }

  /**
   * Circle color: hue encodes direction (red = raises activity, blue = lowers), saturation encodes significance
   * (vivid when highly significant, pale near p≈1, grey when the t-test could not be run). Bounded, so an exact
   * p=0 simply saturates instead of skewing any normalization.
   */
  private contributionColor(contribution: number, pValue: number | null): string {
    if (pValue == null)
      return 'rgb(200, 200, 200)';
    const [r, g, b] = contribution >= 0 ? [198, 40, 40] : [21, 101, 192];
    const sig = Math.min(1, Math.max(0, -Math.log10(Math.max(pValue, 1e-12)) / 3)); // p≤0.001 → full color
    const mix = (c: number): number => Math.round(255 + (c - 255) * sig);
    return `rgb(${mix(r)}, ${mix(g)}, ${mix(b)})`;
  }

  /** Text color for a monomer, taken from its monomer-library background color (falls back to dark grey). */
  private monomerTextColor(monomer: string): string {
    const cached = this._monomerColorCache.get(monomer);
    if (cached != null)
      return cached;
    let color = '#333333';
    try {
      const colors = PeptideUtils.getMonomerLib().getMonomerColors(this._biotype, monomer);
      if (colors?.backgroundcolor)
        color = colors.backgroundcolor;
    } catch {
      // monomer library not ready — keep the neutral fallback
    }
    this._monomerColorCache.set(monomer, color);
    return color;
  }

  /** Compact signed contribution label, e.g. "+1.24" or "-0.4". */
  private formatContribution(v: number): string {
    const s = Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2);
    return `${v >= 0 ? '+' : ''}${s}`;
  }

  /** Truncates text with an ellipsis so it fits within maxWidth on the given canvas context. */
  private fitText(g: CanvasRenderingContext2D, text: string, maxWidth: number): string {
    if (maxWidth <= 0 || g.measureText(text).width <= maxWidth)
      return text;
    let t = text;
    while (t.length > 1 && g.measureText(`${t}…`).width > maxWidth)
      t = t.slice(0, -1);
    return t.length > 1 ? `${t}…` : t;
  }

  /**
   * Shows the invariant-map style activity-distribution tooltip for a basis cell (monomer-at-position vs the rest).
   * @param gridCell - Hovered grid cell.
   * @param x - Tooltip x coordinate.
   * @param y - Tooltip y coordinate.
   * @return - Whether the tooltip was handled.
   */
  private showBasisTooltip(gridCell: DG.GridCell, x: number, y: number): boolean {
    const posName = gridCell.tableColumn?.name;
    if (!gridCell.isTableCell || posName == null || !this._basisPositionSet.has(posName) ||
      gridCell.tableRowIndex == null)
      return false;
    const cand = this._basisByPosition[posName]?.[gridCell.tableRowIndex] ?? null;
    if (cand == null || !cand.monomer || this._scaledActivityCol == null)
      return false;

    const monomerPosition = {monomerOrCluster: cand.monomer, positionOrClusterType: posName};
    highlightMonomerPosition(monomerPosition, this.dataFrame, this.monomerPositionStats);
    // Explain the glyph encoding at the top of the distribution panel.
    const direction = cand.contribution >= 0 ? 'raises' : 'lowers';
    const additionalStats = {
      'Monomer': gridCell.cell.valueString,
      'Circle size': `effect magnitude |Δ| = ${Math.abs(cand.contribution).toFixed(2)}`,
      'Circle color': cand.pValue == null ? `${direction} activity (significance: n/a — too few samples)` :
        `${direction} activity, significance p = ` +
        `${cand.pValue < 0.001 ? cand.pValue.toExponential(1) : cand.pValue.toFixed(3)}`,
    };
    return showTooltip(this.dataFrame, this._scaledActivityCol, [], {
      monomerPosition, x, y, mpStats: this.monomerPositionStats, fromViewer: true, additionalStats,
    });
  }

  /**
   * Runs the additive-model beam search and returns the generated-peptides dataframe.
   * @return - Dataframe with generated peptides and per-peptide statistics.
   */
  generatePeptides(): DG.DataFrame {
    const seqCol = this.dataFrame.getCol(this.sequenceColumnName);
    const scaledActivityCol = scaleActivity(this.dataFrame.getCol(this.activityColumnName), this.activityScaling);
    this._scaledActivityCol = scaledActivityCol; // reused by the basis-cell distribution tooltips
    // Baseline must match the population the statistics were computed over (filtered rows when a filter is active).
    const globalMean = this.dataFrame.filter.anyFalse ?
      DG.Stats.fromColumn(scaledActivityCol, this.dataFrame.filter).avg : scaledActivityCol.stats.avg;

    const positionColumns = this.getPositionColumns();
    const stats: MonomerPositionStats = this.monomerPositionStats;

    const isHigh = this.activityTarget === C.ACTIVITY_TARGET.HIGH;
    // Ordered per-position candidate lists.
    const perPositionCandidates: Candidate[][] = [];
    for (const posCol of positionColumns) {
      const posName = posCol.name;
      const posStats = stats[posName] as PositionStats | undefined;
      if (!posStats)
        continue;

      const monomerEntries: Candidate[] = [];
      let bestFallback: Candidate | null = null; // most-supported monomer, used when nothing passes the filters
      let hasSignificant = false;
      for (const monomer of Object.keys(posStats)) {
        if (monomer === 'general')
          continue;
        const item = posStats[monomer] as StatsItem;
        // Empirical-Bayes shrinkage: damp the marginal effect of monomers with little supporting data so that
        // noisy per-position maxima do not sum into wildly out-of-range predictions (winner's curse).
        const shrinkFactor = item.count / (item.count + this.shrinkageStrength);
        const contribution = (item.mean - globalMean) * shrinkFactor;
        const significant = item.pValue != null && item.pValue <= 0.05;
        hasSignificant ||= significant;
        const candidate: Candidate = {
          position: posName, monomer, contribution, count: item.count, pValue: item.pValue,
          significant, supported: item.count >= this.minSupport,
        };
        if (bestFallback === null || candidate.count > bestFallback.count)
          bestFallback = candidate;
        const passesSupport = item.count >= this.minSupport;
        const passesPValue = this.maxPValue >= 1 || (item.pValue != null && item.pValue <= this.maxPValue);
        if (passesSupport && passesPValue)
          monomerEntries.push(candidate);
      }

      if (monomerEntries.length === 0) {
        // Keep the position fixed to its most common monomer so generated sequences stay aligned and valid.
        if (bestFallback !== null)
          perPositionCandidates.push([bestFallback]);
        continue;
      }

      // When requested, don't vary positions without a significant signal — pin them to the best supported monomer.
      if (this.onlySignificantPositions && !hasSignificant) {
        perPositionCandidates.push([bestFallback ?? monomerEntries[0]]);
        continue;
      }

      // Rank towards the target and keep the top-K.
      monomerEntries.sort((a, b) => isHigh ? b.contribution - a.contribution : a.contribution - b.contribution);
      perPositionCandidates.push(monomerEntries.slice(0, this.candidatesPerPosition));
    }

    const generated = this.beamSearch(perPositionCandidates, isHigh);
    // Position order of the beam-search candidate slots, used to lay out the per-position basis columns.
    const positionOrder = perPositionCandidates.map((candidates) => candidates[0].position);
    return this.buildResultDataFrame(generated, seqCol, globalMean, positionOrder);
  }

  /**
   * Bounded beam search over the per-position candidate monomers, keeping the highest (or lowest) scoring partials.
   * @param perPositionCandidates - Ordered candidate monomers for each position.
   * @param isHigh - Whether higher scaled activity is the objective.
   * @return - Ranked list of complete peptides.
   */
  private beamSearch(perPositionCandidates: Candidate[][], isHigh: boolean): BeamItem[] {
    // Enough breadth to still satisfy peptideCount after novelty filtering.
    const width = Math.max(this.beamWidth, this.peptideCount * (this.excludeExisting ? 2 : 1));
    const better = (a: number, b: number): number => isHigh ? b - a : a - b;

    let beam: BeamItem[] = [{monomers: [], candidates: [], cumContribution: 0}];
    for (const candidates of perPositionCandidates) {
      const expanded: BeamItem[] = [];
      for (const item of beam) {
        for (const candidate of candidates) {
          expanded.push({
            monomers: [...item.monomers, candidate.monomer],
            candidates: [...item.candidates, candidate],
            cumContribution: item.cumContribution + candidate.contribution,
          });
        }
      }
      expanded.sort((a, b) => better(a.cumContribution, b.cumContribution));
      beam = expanded.slice(0, width);
    }
    return beam;
  }

  /**
   * Assembles the generated peptides into a grid-ready dataframe with sequence, prediction and supporting statistics.
   * @param generated - Ranked generated peptides.
   * @param seqCol - Source sequence column (used for notation and novelty checks).
   * @param globalMean - Global mean of the scaled activity.
   * @return - Result dataframe.
   */
  private buildResultDataFrame(generated: BeamItem[], seqCol: DG.Column<string>, globalMean: number,
    positionOrder: string[]): DG.DataFrame {
    const seqHelper = PeptideUtils.getSeqHelper();
    const sh: ISeqHandler = seqHelper.getSeqHandler(seqCol);
    const gapOriginal = GapOriginals[sh.notation as NOTATION] ?? GAP_SYMBOL;
    const joiner = sh.joiner;

    // helm joiner does not add v.20 postfix, so we make sure we remove them
    const existing = new Set<string>(sh.isHelm() ? seqCol.categories.map((s) => s?.endsWith('V2.0') || s?.endsWith('V1.0') ? s.substring(0, s.length - 4) : s) :seqCol.categories);

    const seqStrings: string[] = [];
    const predicted: number[] = [];
    const confidence: number[] = [];
    const minSupportArr: number[] = [];
    const meanPValue: number[] = [];
    const significantPositions: number[] = [];
    const novel: boolean[] = [];
    const basis: string[] = [];
    // Per surviving row: the chosen candidate at each position (aligned to positionOrder).
    const basisRows: Candidate[][] = [];
    for (const item of generated) {
      const seq = joiner(new StringListSeqSplitted(item.monomers, gapOriginal));
      const isNovel = !existing.has(seq);
      if (this.excludeExisting && !isNovel)
        continue;
      if (seqStrings.length >= this.peptideCount)
        break;

      const varied = item.candidates.filter((c) => c.supported || c.significant);
      const consideredCount = Math.max(varied.length, 1);
      let nSignificant = 0;
      let nSupported = 0;
      let minCount = Number.POSITIVE_INFINITY;
      let pValueSum = 0;
      let pValueCount = 0;
      for (const candidate of item.candidates) {
        if (candidate.significant)
          ++nSignificant;
        if (candidate.supported)
          ++nSupported;
        minCount = Math.min(minCount, candidate.count);
        if (candidate.pValue != null) {
          pValueSum += candidate.pValue;
          ++pValueCount;
        }
      }

      const topContribs = [...item.candidates]
        .sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution)).slice(0, 5);
      const basisStr = topContribs
        .map((c) => `${c.position}:${c.monomer}(${c.contribution >= 0 ? '+' : ''}${c.contribution.toFixed(2)})`)
        .join(', ');

      seqStrings.push(seq);
      predicted.push(globalMean + item.cumContribution);
      confidence.push(nSupported === 0 ? 0 : (nSignificant / consideredCount));
      minSupportArr.push(minCount === Number.POSITIVE_INFINITY ? 0 : minCount);
      meanPValue.push(pValueCount === 0 ? DG.FLOAT_NULL : pValueSum / pValueCount);
      significantPositions.push(nSignificant);
      novel.push(isNovel);
      basis.push(basisStr);
      basisRows.push(item.candidates);
    }

    const peptideCol = sh.getNewColumnFromList(GEN_COLUMN_NAMES.PEPTIDE, seqStrings);
    const columns: DG.Column[] = [
      peptideCol,
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, GEN_COLUMN_NAMES.PREDICTED, predicted),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, GEN_COLUMN_NAMES.CONFIDENCE, confidence),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, GEN_COLUMN_NAMES.MIN_SUPPORT, minSupportArr),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, GEN_COLUMN_NAMES.MEAN_P_VALUE, meanPValue),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, GEN_COLUMN_NAMES.SIGNIFICANT_POSITIONS, significantPositions),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, GEN_COLUMN_NAMES.NOVEL, novel),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, GEN_COLUMN_NAMES.BASIS, basis),
    ];

    // One string column per position holding each row's chosen monomer, custom-rendered as monomer + circle + number.
    this._basisByPosition = {};
    this._basisPositionNames = [];
    this._rowMaxAbsContribution = basisRows.map((cands) =>
      cands.reduce((mx, c) => Math.max(mx, Math.abs(c.contribution)), 0));
    const usedNames = new Set(columns.map((col) => col.name));
    for (let j = 0; j < positionOrder.length; ++j) {
      let colName = positionOrder[j];
      while (usedNames.has(colName)) // avoid clashing with the statistics columns
        colName = `${colName} `;
      usedNames.add(colName);
      const cellCandidates = basisRows.map((cands) => cands[j] ?? null);
      columns.push(DG.Column.fromStrings(colName, cellCandidates.map((c) => c?.monomer ?? '')));
      this._basisByPosition[colName] = cellCandidates;
      this._basisPositionNames.push(colName);
    }
    this._basisPositionSet = new Set(this._basisPositionNames);

    const resultDf = DG.DataFrame.fromColumns(columns);
    resultDf.name = 'Generated Peptides';
    return resultDf;
  }
}
