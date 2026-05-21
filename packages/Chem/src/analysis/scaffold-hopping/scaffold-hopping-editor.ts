import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getRdKitModule, RDKIT_COMMON_RENDER_OPTS} from '../../utils/chem-common-rdkit';
import {extractAtomPositionsFromSvg, findNearestAtom, computeSelectedBonds}
  from '../../utils/chem-atom-picker-utils';
import * as chemSearches from '../../chem-searches';

/** RGB triple in [0, 1] used by RDKit's `get_svg_with_highlights` for the
 *  user-marked replaceable region. Orange — visually distinct from the default
 *  RDKit query/match colours. */
const REPLACEABLE_COLOR: [number, number, number] = [1.0, 0.6, 0.2];

/** SVG canvas size for the reference preview. Square; RDKit will fit the
 *  molecule into the box and the picker uses these as click coordinates. */
const PREVIEW_SIZE = 320;

/** Four named threshold presets. **Category labels follow the qualitative
 *  Sun-Tawa-Wallqvist 2012 hop taxonomy** ("Classification of scaffold
 *  hopping approaches", DDT 17:310), which describes the structural-change
 *  regimes by name but does NOT prescribe numeric cutoffs. **The numeric
 *  thresholds below are Datagrok defaults**, calibrated on the BCR-ABL
 *  smoke set and tuned for usable defaults across the four regimes; only
 *  the Hard preset's `mcsRatioMax = 0.4` matches a published number
 *  (Maeda 2024's atom-ratio cap for scaffold hops). Users are expected to
 *  edit the numeric inputs in the dialog when their dataset's distribution
 *  warrants it.
 *
 *  - **Local** (R-group / bioisostere refinement, marked-atoms driven):
 *    user marks a sub-region of the reference and the run looks for
 *    candidates that share the unmarked region but swap the marked part.
 *    Tc 0.5-0.95 keeps the candidate inside the reference's chemotype;
 *    mcsRatio ≤ 0.98 admits single-atom bioisosteric swaps (O↔N in
 *    bicyclic rings, Cl↔F, ring-N position changes) that 0.95 silently
 *    dropped on drug-sized molecules (a 30-atom molecule with one swap
 *    gives ratio 29/30 = 0.97). The marked-atoms criterion remains the
 *    load-bearing filter; mcsRatio is the backstop against zero-swap
 *    near-duplicates. Datagrok defaults; not held-out-validated.
 *  - **Easy** (substituent / heterocycle swap, Sun-Tawa-Wallqvist class 1°):
 *    candidates share most of the scaffold; hops are decorations or single-
 *    ring isosteres. Tc 0.4-0.7 is the typical range; mcsRatio ≤ 0.75 admits
 *    single-heteroaryl-ring swaps on medium drug-sized molecules
 *    (a 40-atom molecule with a 10-atom ring isostere gives ratio ≈ 0.75).
 *    Both numbers are Datagrok defaults.
 *  - **Middle** (ring open/close, residue swap, class 2°): meaningful but
 *    bounded scaffold change. Default, balanced for typical SAR-table
 *    workflows. Datagrok defaults.
 *  - **Hard / Maeda-faithful** (large topological change, class 3°/4°):
 *    paper-faithful Maeda 2024 atom-ratio cap (0.4). The Tc lower bound
 *    (0.05) is a Datagrok default informed by Schneider 1999, Vogt et al.
 *    2010, and Iktos/Pinel 2023, which all show genuine large-topology hops
 *    sit in the 0.05-0.30 ECFP4 range — higher cutoffs systematically
 *    discard exactly the candidates Maeda's classifier is designed to find. */
type PresetKey = 'Local' | 'Easy' | 'Middle' | 'Hard';
// Per-preset defaults. `activityWeight` is the proximity-multiplier
// strength in [0, 1]: Local hops share a SARM with the reference so MMP
// predictions are well-calibrated and a stronger weight (0.30) earns its
// keep; non-Local presets cross chemotypes, where both MMP and kNN
// predictions are intrinsically noisier, so the activity nudge is kept
// softer (0.20 Easy → 0.15 Middle/Hard) to avoid noisy predictions
// rewriting the structural ranking too aggressively. Proximity stays on
// by default everywhere — turning it off entirely would hide the feature
// from one-shot users and waste the MMPA/kNN compute cost.
const PRESETS: Record<PresetKey, {tcMin: number; tcMax: number; mcsRatioMax: number;
  activityWeight: number}> = {
  // Local mcsRatioMax = 0.98 (was 0.95): on drug-sized molecules (~30-50
  // heavy atoms) a single-atom swap such as O↔N in a bicyclic ring,
  // Cl↔F on a substituent, or C↔N in a pyridine→pyrimidine bioisostere
  // gives mcsRatio = 29/30 to 49/50 ≈ 0.967-0.980. At 0.95 those swaps
  // sat just above the cap and were never flagged as hops — exactly the
  // bioisosteric class Local mode is supposed to surface. 0.98 admits
  // single-atom heteroatom swaps; identical molecules (ratio = 1.0)
  // still get filtered out. The marked-atoms criterion remains the
  // load-bearing filter when atoms are marked. NOT held-out-validated.
    'Local': {tcMin: 0.5, tcMax: 0.95, mcsRatioMax: 0.98, activityWeight: 0.3},
    // Easy mcsRatioMax = 0.75 (was 0.70): single-ring isostere swaps on
    // medium-sized drugs (40-atom molecule with a 10-atom heteroaryl ring
    // replaced by an isosteric ring of similar size) give mcsRatio ≈ 0.75
    // — exactly at the old cap. 0.70 systematically dropped them. 0.75
    // still requires ≥25% atom-level structural change, which keeps the
    // class-1° "substituent / heterocycle swap" semantics intact. NOT
    // held-out-validated; rationale matches the Local 0.95→0.98 change.
    'Easy': {tcMin: 0.4, tcMax: 0.7, mcsRatioMax: 0.75, activityWeight: 0.2},
    'Middle': {tcMin: 0.2, tcMax: 0.5, mcsRatioMax: 0.5, activityWeight: 0.15},
    'Hard': {tcMin: 0.05, tcMax: 0.3, mcsRatioMax: 0.4, activityWeight: 0.15},
  };
/** Display labels for the preset dropdown — keep the taxonomy hint visible
 *  so the user doesn't have to read the tooltip to understand the choice.
 *  Local appears first as the most directed (region-specific) mode; the
 *  Easy/Middle/Hard sequence after that is in increasing structural-
 *  distance order so the user can scan top-to-bottom by "how far do I
 *  want to hop?" */
const PRESET_LABELS: Record<PresetKey, string> = {
  'Local': 'Local — change only the marked region (requires atom marks)',
  'Easy': 'Easy — substituent / heterocycle swap',
  'Middle': 'Middle — ring open/close, residue swap',
  'Hard': 'Hard — large topological change (atom-ratio ≤ 0.4)',
};
const DEFAULT_PRESET: PresetKey = 'Middle';

/** localStorage key for cross-session preset persistence. The chosen
 *  preset survives full page reloads — first-run users see Middle; users
 *  who've picked Hard previously see Hard on every subsequent visit. */
const PRESET_STORAGE_KEY = 'datagrok.chem.scaffoldHopping.lastPreset';

/** Reads the last-selected preset from localStorage, falling back to the
 *  default when nothing is stored or the stored value is unrecognised. */
function loadLastPreset(): PresetKey {
  try {
    const stored = (typeof localStorage !== 'undefined') ?
      localStorage.getItem(PRESET_STORAGE_KEY) : null;
    if (stored && stored in PRESETS) return stored as PresetKey;
  } catch {/* localStorage disabled or full — fall through */}
  return DEFAULT_PRESET;
}

/** Persists the user's preset choice to localStorage. Non-fatal on failure
 *  (private-mode browsers, quota exhausted, etc.) — the module-level
 *  in-session cache still keeps the choice alive for the current page. */
function persistLastPreset(key: PresetKey): void {
  try {
    if (typeof localStorage !== 'undefined')
      localStorage.setItem(PRESET_STORAGE_KEY, key);
  } catch {/* ignore */}
}

/** Module-level cache of the user's last-selected preset. Persists for the
 *  page lifetime so the dialog remembers the choice across open/close
 *  cycles within a session. Seeded from `loadLastPreset()` on module load
 *  so the FIRST dialog open in a new session also gets the persisted
 *  value, not just subsequent opens. Updated by `_applyPreset` whenever
 *  the user picks something. */
let _lastPresetKey: PresetKey = loadLastPreset();

/** Common activity-column name patterns the auto-detector recognises.
 *  Case-insensitive substring match against `column.name`. First match
 *  wins. The list intentionally covers both p-scaled and raw scales —
 *  the proximity computation is scale-agnostic (sigma is derived from
 *  the column's own std-dev). */
const ACTIVITY_COLUMN_PATTERNS = [
  'pic50', 'pki', 'pec50', 'pkd', 'pactivity',
  'ic50', 'ki', 'ec50', 'kd',
  'activity', 'potency', 'affinity',
];

/** Module-level cache of replaceable-region atom marks, keyed by reference
 *  SMILES. Persists for the page lifetime — across dialog open/close cycles
 *  and across different tables that share the same reference molecule. The
 *  user only has to mark a region once per molecule; subsequent invocations
 *  restore the marks automatically. Cleared by the user via the Clear-
 *  selection button (which also clears the cache entry for that SMILES). */
const _atomMarkCache = new Map<string, number[]>();

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
  /** Activity-aware re-rank — empty when disabled. When set, names a
   *  numeric column in `table`; each candidate's composite score is
   *  multiplied by (1 - activityWeight) + activityWeight × proximity. */
  activityColumnName: string;
  /** Blend strength of the activity proximity multiplier, in [0, 1]. */
  activityWeight: number;
  /** When true, the run writes an extra `Scaffold Hop Reason` string
   *  column with the per-row pass/fail breakdown. Off by default to keep
   *  the result table compact; chemists who want to inspect why a given
   *  row was or wasn't flagged can toggle it on. */
  showReason: boolean;
  /** When true, the per-row Replacement / Replaced Region columns are
   *  computed via R-group decomposition instead of ErG pharmacophore
   *  matching. Set automatically when the user picks the Local preset. */
  useRGroupReplacement: boolean;
  /** When true, run MMP-rule-based activity prediction on
   *  `activityColumnName` (same machinery the MMP viewer's Generations tab
   *  uses). Adds `<actName> (MMP Pred / Support / Stdev)` columns and lets
   *  the proximity loop fill in masked / unknown activities. */
  imputeActivity: boolean;
  /** Minimum supporting pairs per MMP rule (default 10). */
  mmpSupportFloor: number;
  /** Also predict activity for rows that already have a measured value
   *  (validation mode). Default false. */
  predictKnown: boolean;
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
  /** Either a `choice<string>` (when the table has a name-like column —
   *  shows "name (row N)" labels, sorted) or an `int` (fallback when no
   *  name column exists). Read the row index via `_getCurrentRefIdx`,
   *  which handles both shapes. */
  refRowInput!: DG.InputBase<any>;
  /** When the choice variant is in use, maps each label string back to
   *  its 0-based row index. `null` when the int fallback is in use. */
  private _refIdxLookup: Map<string, number> | null = null;

  presetInput = ui.input.choice<string>('Preset', {
    value: PRESET_LABELS[_lastPresetKey],
    items: (Object.keys(PRESETS) as PresetKey[]).map((k) => PRESET_LABELS[k]),
    onValueChanged: () => this._applyPreset(),
    tooltipText: 'Picks Tc range and MCS atom-ratio cap to match the kind of ' +
      'scaffold change you\'re looking for. Pick Easy for substituent / ' +
      'heterocycle swaps, Middle for ring open/close or residue swaps, Hard ' +
      'for large topological changes (strict atom-ratio ≤ 0.4). The three ' +
      'numeric inputs below stay editable — the preset just seeds them.',
  });

  tanimotoMinInput = ui.input.float('ECFP4 Tc min', {value: PRESETS[_lastPresetKey].tcMin,
    tooltipText: 'Pre-filter shortlist lower bound. Drops rows with negligible ' +
      'similarity before pharmacophore / MCS scoring. Lower this if your ' +
      'reference and the genuine hops are expected to share little topology — ' +
      'Schneider 1999 and Iktos/Pinel 2023 both report real hops at Tc < 0.2.'});
  tanimotoMaxInput = ui.input.float('ECFP4 Tc max', {value: PRESETS[_lastPresetKey].tcMax,
    tooltipText: 'Pre-filter shortlist upper bound. Drops near-duplicates of ' +
      'the reference (Tc above this are usually trivial analogues, not hops).'});
  mcsMaxInput = ui.input.float('MCS atom ratio max', {value: PRESETS[_lastPresetKey].mcsRatioMax,
    tooltipText: 'Scaffold-hop atom-ratio threshold: a candidate is flagged ' +
      'as a hop iff ratio_atom = atoms(MCS) / atoms(reference) is ≤ this ' +
      'threshold. Hard preset uses 0.4 (strict — at least 60% of the ' +
      'reference must change); Middle and Easy presets relax it.'});
  useTcInFlagInput = ui.input.bool('Tc window in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires the ' +
      'candidate\'s ECFP4 Tc to be inside [Tc min, Tc max]. When unchecked, ' +
      'the Tc range is only used as a pre-filter shortlist for performance — ' +
      'the flag itself ignores it. Uncheck this and "CATS Sim in flag" ' +
      'together to fall back to the strict atom-ratio definition (≤ 0.4 only).'});
  useCatsInFlagInput = ui.input.bool('CATS Sim in flag', {value: true,
    tooltipText: 'When checked, the Scaffold Hop flag also requires CATS Sim ' +
      '≥ CATS Sim min. When unchecked, CATS similarity is computed and shown ' +
      'but does not affect the flag — leaves only the strict atom-ratio ' +
      'criterion (≤ 0.4) driving the flag.'});

  showReasonInput = ui.input.bool('Show reason column', {value: false,
    tooltipText: 'When checked, the result table gets an extra ' +
      '"Scaffold Hop Reason" string column with a per-row pass/fail breakdown ' +
      'of every flag condition (e.g. "Tc 0.18 ✓ · CATS 0.81 ✓ · MCS 0.52 ✗"). ' +
      'Useful when you want to know exactly WHY a candidate was or wasn\'t ' +
      'flagged. Off by default to keep the result table compact.'});

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

  /** Activity-column dropdown — choices are the table's numeric columns
   *  matching the activity-name patterns, plus "(none)" to disable
   *  activity weighting. Built dynamically in `onTableInputChanged` so
   *  it tracks the current table. */
  activityColInput!: DG.InputBase<string | null>;
  activityColInputRoot: HTMLElement = ui.div();
  activityWeightInput = ui.input.float('Activity weight', {value: PRESETS[_lastPresetKey].activityWeight,
    tooltipText: 'Strength of the activity-proximity multiplier on the ' +
      'composite score, in [0, 1]. 0 = activity is computed and shown in a ' +
      'separate column but does NOT shift the ranking — recommended for ' +
      'prospective discovery where activities aren\'t yet known. The preset ' +
      'defaults nudge activity-similar candidates toward the top without ' +
      'losing the structural signal: 0.30 for Local (well-calibrated within ' +
      'a shared scaffold series), 0.15-0.20 for non-Local presets where ' +
      'cross-chemotype predictions are noisier. The proximity itself is ' +
      'exp(-|delta| / sigma) where sigma = the activity column\'s standard ' +
      'deviation. EXPECTS LOG-SCALED ACTIVITY (pIC50 / pKi / pEC50). Raw ' +
      'IC50/Ki/Kd columns are log-normal-distributed; the wide range blows ' +
      'sigma out so exp(-|delta|/sigma) collapses to ~1 for every pair and ' +
      'the re-rank does nothing. Pre-compute -log10(activity) into a new ' +
      'column and pick that instead if your data is raw.'});

  /** Optional MMP-rule-based prediction of missing activities. Same
   *  machinery the MMP viewer's Generations tab uses — `MMPA.init` to
   *  build the rule set, then aggregate `meanDiff` over anchor pairs to
   *  predict each row's value. Falls back to kNN-on-Morgan-FP for rows
   *  MMP can't cover. Adds three columns to the table and lets the
   *  proximity loop fill in masked rows. Off by default — running MMP
   *  analysis on every scaffold-hop run is expensive on large tables. */
  imputeActivityInput = ui.input.bool('Predict missing activity (MMP)', {value: false,
    tooltipText: 'When checked, runs MMP-rule-based activity prediction on ' +
      'the selected activity column with a Morgan-FP kNN fallback. Adds ' +
      'three Scaffold Hop Predicted Activity columns and fills masked / ' +
      'unknown activities for the proximity-weighted ranking.'});

  /** Minimum supporting pairs per MMP rule. mmpdb (Dalke 2018) prefers
   *  >=10; relaxing to 5 is its documented fallback. Higher = fewer but
   *  more reliable rules (Kramer 2014's significance gate tightens at
   *  n>=10); lower = more coverage but noisier per-anchor meanDiffs. */
  mmpSupportFloorInput = ui.input.int('MMP support floor', {value: 10,
    tooltipText: 'Minimum number of supporting molecule pairs an MMP rule ' +
      'must have before it can produce a prediction. Default 10 matches ' +
      'mmpdb\'s preferred threshold (Dalke 2018). Drop to 5 for sparse ' +
      'datasets — more MMP coverage at the cost of noisier per-rule ' +
      'meanDiffs. Below 5 the predictions stop being meaningful.'});

  /** Whether MMP predictions are also computed for rows with measured
   *  activity. Default off in production — emitting a prediction for a
   *  row whose activity is already known introduces mild self-inclusion
   *  bias (the rule's meanDiff includes the target's own pairs). Turn on
   *  to compare prediction vs measured side by side for validation. */
  predictKnownInput = ui.input.bool('Predict measured rows too (validation)', {value: false,
    tooltipText: 'When checked, MMP also predicts activity for rows that ' +
      'have a measured value, so you can compare prediction vs measured ' +
      'side by side. Off by default — including the target\'s own pairs ' +
      'in the rule meanDiff introduces a mild self-inclusion bias that\'s ' +
      'fine for sanity-checking but inappropriate for production runs.'});

  /** Live row-count preview shown beneath the Tc inputs: "~120 of 4,217 ' +
   *  rows will pass pre-filter". Updated on every Tc-input or reference
   *  change; uses a 200-row Tanimoto sample (cached in `_sampleTc`) so
   *  the update is sub-second even on 100k-row tables. */
  rowCountPreviewLabel: HTMLElement = ui.divText('', {style: {
    fontSize: '11px', color: 'var(--grey-6)', marginTop: '4px',
    fontStyle: 'italic',
  }});
  /** Cached Tc values for a random sample of rows vs. the current
   *  reference. Null until the first compute completes. Recomputed
   *  whenever the reference row or molecule column changes. */
  private _sampleTc: Float32Array | null = null;
  /** Row indices that were sampled — same length as `_sampleTc`. Used
   *  only to keep the sample-vs-total ratio honest in the preview. */
  private _sampleRowIdxs: number[] = [];

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
  /** Permanent modifier-key hint, always visible below the preview so the
   *  Ctrl-hover-paint affordance is discoverable on first use without
   *  having to mouse over a tooltip. */
  modifierHintLabel: HTMLElement = ui.divText('', {style: {
    fontSize: '11px', color: 'var(--grey-5)', marginTop: '6px',
    maxWidth: `${PREVIEW_SIZE}px`, lineHeight: '1.4',
  }});
  selectionCaption: HTMLElement = ui.divText('', {style: {
    fontSize: '11px', color: 'var(--grey-6)', marginTop: '4px', fontWeight: '500',
    maxWidth: `${PREVIEW_SIZE}px`,
  }});
  clearSelectionBtn = ui.button('Clear selection', () => this._clearSelection(),
    'Remove all marked atoms — the run will return any global scaffold hop');

  /** Live "what will run" summary at the top of the right panel. Reads
   *  the current table / reference / preset / threshold inputs and
   *  rewrites itself on every change — gives the user a single sanity-
   *  check sentence before they click OK ("am I scoring the right
   *  reference against the right candidates with the right thresholds?"). */
  summaryLabel: HTMLElement = ui.divText('', {style: {
    fontSize: '12px', color: 'var(--grey-7)', fontStyle: 'italic',
    padding: '6px 0', borderBottom: '1px solid var(--grey-2)', marginBottom: '8px',
  }});

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
    this.activityWeightInput.addValidator(this._rangeValidator);
    // Wire the live summary to every threshold input — preset changes
    // also trigger this via `_applyPreset` which writes to all three.
    // The Tc inputs also drive the live row-count preview, so they're
    // double-subscribed.
    this.tanimotoMinInput.onChanged.subscribe(() => {
      this._updateSummary(); this._updateRowCountPreview();
    });
    this.tanimotoMaxInput.onChanged.subscribe(() => {
      this._updateSummary(); this._updateRowCountPreview();
    });
    this.mcsMaxInput.onChanged.subscribe(() => this._updateSummary());
    this.catsSimInput.onChanged.subscribe(() => this._updateSummary());
    this.activityWeightInput.onChanged.subscribe(() => this._updateSummary());
    this.modifierHintLabel.innerHTML =
      '<span style="font-weight:500">Mark replaceable region:</span> ' +
      'Ctrl-hover to paint • Ctrl+Shift-hover to erase • ' +
      '<span style="color:var(--grey-4)">leave empty for any global hop</span>';
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
      onValueChanged: () => {this._refreshReferencePreview(); this._updateSummary();},
    });
    ui.empty(this.colInputRoot);
    Array.from(newColInput.root.children).forEach((c) => this.colInputRoot.append(c));
    this.colInput = newColInput;

    const defaultIdx = table.currentRowIdx >= 0 ? table.currentRowIdx : 0;
    const nameCol = this._findNameLikeColumn(table, firstMol);
    if (nameCol) {
      // Choice variant: "name (row N)" labels, sorted alphabetically by
      // name. Order of magnitude more usable than typing row indices for
      // any table that has a name column (which is most chemists' SAR
      // tables in practice).
      const labelByIdx: string[] = [];
      const lookup = new Map<string, number>();
      for (let i = 0; i < table.rowCount; i++) {
        const raw = nameCol.get(i);
        const display = (typeof raw === 'string' && raw.length > 0) ? raw : '(unnamed)';
        const label = `${display} (row ${i})`;
        labelByIdx.push(label);
        lookup.set(label, i);
      }
      const sortedLabels = labelByIdx.slice().sort((a, b) =>
        a.localeCompare(b, undefined, {numeric: true, sensitivity: 'base'}));
      this._refIdxLookup = lookup;
      this.refRowInput = ui.input.choice<string>('Reference row', {
        value: labelByIdx[defaultIdx] ?? sortedLabels[0],
        items: sortedLabels,
        onValueChanged: () => {this._refreshReferencePreview(); this._updateSummary();},
      });
    } else {
      // Fallback: integer input with the original range validator. Only
      // used when the table has no name-like column at all.
      this._refIdxLookup = null;
      this.refRowInput = ui.input.int('Reference row', {
        value: defaultIdx,
        onValueChanged: () => {this._refreshReferencePreview(); this._updateSummary();},
      });
      this.refRowInput.addValidator((s: string) => {
        const v = parseInt(s, 10);
        if (Number.isNaN(v)) return 'Integer expected';
        if (v < 0 || v >= table.rowCount) return `Must be in [0, ${table.rowCount - 1}]`;
        return null;
      });
    }

    // Activity-column picker. Filter to numeric columns whose name matches
    // a known activity pattern. "(none)" is always the first item so the
    // user can disable activity weighting explicitly. The auto-detected
    // best match becomes the default selection.
    this._rebuildActivityColumnInput(table);

    this._refreshReferencePreview();
    this._updateSummary();
  }

  /** Rebuilds the activity-column dropdown for the current table. The
   *  choices are: "(none)" + every numeric column whose name matches one
   *  of `ACTIVITY_COLUMN_PATTERNS`. When at least one activity column is
   *  detected, the first match wins as the default; otherwise the input
   *  defaults to "(none)" and activity weighting stays off until the
   *  user picks one. */
  private _rebuildActivityColumnInput(table: DG.DataFrame): void {
    const numericCols = table.columns.toList().filter((c) =>
      c.type === DG.COLUMN_TYPE.FLOAT || c.type === DG.COLUMN_TYPE.INT ||
      c.type === DG.COLUMN_TYPE.BIG_INT);
    const detected: string[] = [];
    const other: string[] = [];
    for (const c of numericCols) {
      const lower = c.name.toLowerCase();
      if (ACTIVITY_COLUMN_PATTERNS.some((p) => lower.includes(p)))
        detected.push(c.name);
      else
        other.push(c.name);
    }
    const items = ['(none)', ...detected, ...other];
    const defaultValue = detected.length > 0 ? detected[0] : '(none)';
    const newInput = ui.input.choice<string>('Activity column', {
      value: defaultValue,
      items,
      onValueChanged: () => this._updateSummary(),
      tooltipText: 'Numeric column carrying the activity readout (pIC50, ' +
        'pKi, IC50, etc.). When set, each candidate\'s composite score is ' +
        'multiplied by an activity-proximity factor so hits near the ' +
        'reference\'s potency rank higher. The factor itself is exposed as a ' +
        'separate column ("Scaffold Hop Activity Proximity") for transparency. ' +
        'Auto-detected from column names matching pIC50 / pKi / pEC50 / IC50 / ' +
        'Ki / EC50 / Kd / activity / potency / affinity (case-insensitive).',
    });
    ui.empty(this.activityColInputRoot);
    Array.from(newInput.root.children).forEach((c) => this.activityColInputRoot.append(c));
    this.activityColInput = newInput;
  }

  /** Resolves the current reference row index regardless of which input
   *  variant is active. The int variant returns `value` directly; the
   *  choice variant looks up the label in `_refIdxLookup`. Returns `null`
   *  when no valid selection exists. */
  private _getCurrentRefIdx(): number | null {
    const v = this.refRowInput?.value;
    if (typeof v === 'number') return Number.isFinite(v) ? v : null;
    if (typeof v === 'string' && this._refIdxLookup) {
      const idx = this._refIdxLookup.get(v);
      return idx === undefined ? null : idx;
    }
    return null;
  }

  /** Returns the first name-like column in `table` other than `molCol`,
   *  or `null` if none exists. Same column-name preference order as
   *  `_referenceDisplayName` so the picker and summary agree on which
   *  column drives the label. */
  private _findNameLikeColumn(
    table: DG.DataFrame, molCol: DG.Column | null | undefined,
  ): DG.Column | null {
    const candidates = ['name', 'Name', 'title', 'Title', 'label', 'Label',
      'id', 'ID', 'Id', 'compound', 'Compound'];
    for (const n of candidates) {
      const c = table.col(n);
      if (c && c !== molCol) return c;
    }
    return null;
  }

  private _refreshReferencePreview() {
    const table = this.tableInput.value;
    const col = this.colInput?.value;
    const rowIdx = this._getCurrentRefIdx();
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
    // Restore previously-marked atoms for this molecule (if any). The cache
    // is keyed by SMILES string so the same reference molecule gets the same
    // marks across dialog open/close cycles, table changes, row reordering,
    // etc. — the user only has to mark a region once per molecule.
    const cached = _atomMarkCache.get(smiles);
    if (cached && cached.length > 0)
      this.selectedAtoms = new Set(cached);
    this.referenceCaption.textContent = `Row ${rowIdx} • ${smiles}`;
    this._renderInteractivePicker(smiles);

    // Kick off (or skip) the sample-Tc computation that drives the live
    // row-count preview. Reference changed → cached sample is stale, so
    // invalidate it before the new compute starts. The compute is async;
    // the preview shows "computing…" until it completes, then re-renders
    // with the actual count.
    this._sampleTc = null;
    this._sampleRowIdxs = [];
    this._updateRowCountPreview();
    void this._refreshSampleTc(table, col, smiles);
  }

  /** Asynchronously computes Tanimoto similarities for a random 200-row
   *  sample of the table vs. the current reference SMILES, then refreshes
   *  the row-count preview. Sample-based to keep the dialog responsive
   *  on very large tables — for tables under the sample size we just
   *  compute against everything. Concurrent updates (user switching
   *  references rapidly) are tolerated via the
   *  `expectedRefSmiles === this._currentRefSmiles` guard before write. */
  private async _refreshSampleTc(
    table: DG.DataFrame, col: DG.Column, expectedRefSmiles: string,
  ): Promise<void> {
    const SAMPLE_SIZE = 200;
    const total = table.rowCount;
    if (total <= 1) return;
    const indices: number[] = [];
    if (total - 1 <= SAMPLE_SIZE)
      for (let i = 0; i < total; i++) indices.push(i);
    else {
      // Reservoir-style uniform sample without replacement. Fast and
      // deterministic enough for a UI preview — proper PRNG seeding is
      // not warranted for an approximate count.
      const taken = new Set<number>();
      while (taken.size < SAMPLE_SIZE) taken.add(Math.floor(Math.random() * total));
      for (const i of taken) indices.push(i);
    }
    const sampleCol = DG.Column.fromStrings('__shp_sample__',
      indices.map((i) => col.get(i) ?? ''));
    sampleCol.semType = DG.SEMTYPE.MOLECULE;
    let tcCol: DG.Column | null = null;
    try {
      tcCol = await chemSearches.chemGetSimilarities(sampleCol, expectedRefSmiles);
    } catch {/* sample compute failed — preview stays "computing…" */}
    if (this._currentRefSmiles !== expectedRefSmiles) return; // user switched ref mid-flight
    if (!tcCol) return;
    const arr = new Float32Array(indices.length);
    for (let i = 0; i < indices.length; i++) arr[i] = tcCol.get(i) ?? NaN;
    this._sampleTc = arr;
    this._sampleRowIdxs = indices;
    this._updateRowCountPreview();
  }

  /** Rewrites the row-count preview label from the cached `_sampleTc`
   *  array plus the current Tc-min / Tc-max input values. Estimate when
   *  the sample is a proper subset; exact when the sample is the whole
   *  table. Shows "computing…" until the first sample arrives. */
  private _updateRowCountPreview(): void {
    const table = this.tableInput.value;
    if (!table || !this._currentRefSmiles) {
      this.rowCountPreviewLabel.textContent = '';
      return;
    }
    const total = table.rowCount;
    if (!this._sampleTc) {
      this.rowCountPreviewLabel.textContent = 'Estimating pre-filter survivors…';
      return;
    }
    const tcMin = this.tanimotoMinInput.value ?? 0;
    const tcMax = this.tanimotoMaxInput.value ?? 1;
    let passed = 0;
    for (let i = 0; i < this._sampleTc.length; i++) {
      const v = this._sampleTc[i];
      if (Number.isNaN(v)) continue;
      if (v >= tcMin && v <= tcMax) passed++;
    }
    const candidateCount = Math.max(0, total - 1);
    if (this._sampleTc.length >= total) {
      // We sampled the whole table — exact count.
      this.rowCountPreviewLabel.textContent =
        `${passed.toLocaleString()} of ${candidateCount.toLocaleString()} rows will pass pre-filter.`;
    } else {
      // Scale sample fraction to the full table.
      const fraction = passed / this._sampleTc.length;
      const estimate = Math.round(fraction * candidateCount);
      this.rowCountPreviewLabel.textContent =
        `~${estimate.toLocaleString()} of ${candidateCount.toLocaleString()} rows will pass pre-filter ` +
        `(sample of ${this._sampleTc.length}).`;
    }
  }

  /** Writes the current `selectedAtoms` to the module-level cache, keyed by
   *  the current reference SMILES. Empty selections are deleted from the
   *  cache so an explicit Clear-selection drops the cache entry too. */
  private _saveMarksToCache() {
    if (!this._currentRefSmiles) return;
    if (this.selectedAtoms.size === 0)
      _atomMarkCache.delete(this._currentRefSmiles);
    else
      _atomMarkCache.set(this._currentRefSmiles, [...this.selectedAtoms]);
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
    if (isErase)
      changed = this.selectedAtoms.delete(nearest);
    else if (isPaint) {
      if (!this.selectedAtoms.has(nearest)) {
        this.selectedAtoms.add(nearest);
        changed = true;
      }
    } else {
      if (this.selectedAtoms.has(nearest)) this.selectedAtoms.delete(nearest);
      else this.selectedAtoms.add(nearest);
      changed = true;
    }
    if (changed) {
      this._saveMarksToCache();
      this._reRender(smiles);
    }
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
    if (mode === 'erase')
      changed = this.selectedAtoms.delete(nearest);
    else {
      if (!this.selectedAtoms.has(nearest)) {
        this.selectedAtoms.add(nearest);
        changed = true;
      }
    }
    if (changed) {
      this._saveMarksToCache();
      this._reRender(smiles);
    }
  }

  private _resetPaintTracking() {
    this._lastPaintedAtom = null;
    this._lastPaintMode = null;
  }

  /** Applies the currently-selected preset's Tc/mcsRatio values to the three
   *  numeric inputs and remembers the choice in the module-level cache so
   *  the next dialog open shows the same preset. The user can still tweak
   *  the inputs manually after the preset is applied — we don't lock them. */
  private _applyPreset() {
    const label = this.presetInput.value;
    const key = (Object.keys(PRESETS) as PresetKey[]).find((k) => PRESET_LABELS[k] === label);
    if (!key) return;
    _lastPresetKey = key;
    const p = PRESETS[key];
    this.tanimotoMinInput.value = p.tcMin;
    this.tanimotoMaxInput.value = p.tcMax;
    this.mcsMaxInput.value = p.mcsRatioMax;
    this.activityWeightInput.value = p.activityWeight;
    // Local preset has three additional behaviours: the Tc-window and
    // CATS-Sim flag toggles get switched off (in Local mode the atom mask
    // is the load-bearing filter — Tc is shortlist only, and a 1–2 atom
    // swap barely shifts CATS), CATS Sim min drops to 0.3 (permissive
    // floor so the score column doesn't push interesting local hops
    // down when CATS dips slightly), and the Reason column is auto-
    // enabled (verifying which atoms the algorithm actually changed is
    // the load-bearing diagnostic for local hops). These are seeds, not
    // locks — the user can still flip any of them after the preset
    // applies.
    if (key === 'Local') {
      this.useTcInFlagInput.value = false;
      this.useCatsInFlagInput.value = false;
      this.catsSimInput.value = 0.3;
      this.showReasonInput.value = true;
    }
    persistLastPreset(key);
    this._updateSummary();
    this._updateRowCountPreview();
  }

  /** Rewrites the live "what will run" summary at the top of the right
   *  panel from the current input values. Called from every input's
   *  onValueChanged. Format: `Will score N candidates against <ref> using
   *  <Preset> preset (Tc <min>–<max>, atom-ratio ≤ <mcs>).` Falls back to
   *  empty when the table / column / refRow aren't valid yet. */
  private _updateSummary(): void {
    const table = this.tableInput.value;
    const col = this.colInput?.value;
    const refIdx = this._getCurrentRefIdx();
    if (!table || !col || refIdx == null || refIdx < 0 || refIdx >= table.rowCount) {
      this.summaryLabel.textContent = '';
      return;
    }
    const candidateCount = Math.max(0, table.rowCount - 1).toLocaleString();
    const refName = this._referenceDisplayName(table, refIdx, col);
    const presetLabel = this.presetInput.value;
    const presetKey = (Object.keys(PRESETS) as PresetKey[])
      .find((k) => PRESET_LABELS[k] === presetLabel) ?? _lastPresetKey;

    // Local mode is meaningless without marked atoms — the flag falls
    // back to "any near-neighbor that's not identical" and admits far
    // too many rows on typical SAR tables. Replace the normal summary
    // with a warning until marks are picked. The warning is soft (not a
    // hard block on OK) because some workflows might genuinely want the
    // permissive interpretation, but the text is orange-tinted so the
    // missing-marks case is hard to miss.
    if (presetKey === 'Local' && this.selectedAtoms.size === 0) {
      this.summaryLabel.textContent =
        '⚠ Local mode requires marked atoms — Ctrl-hover atoms on the reference ' +
        'molecule to mark a region. Without marks, every near-neighbor will be ' +
        'flagged as a hop.';
      this.summaryLabel.style.color = 'var(--orange-2, #c87325)';
      return;
    }
    // Normal (non-warning) styling restored for every other state.
    this.summaryLabel.style.color = 'var(--grey-7)';
    const tcMin = this.tanimotoMinInput.value;
    const tcMax = this.tanimotoMaxInput.value;
    const mcs = this.mcsMaxInput.value;
    const cats = this.catsSimInput.value;
    const fmt = (v: number | null | undefined) =>
      typeof v === 'number' && Number.isFinite(v) ?
        (Number.isInteger(v) ? `${v}` : v.toFixed(2)) : '?';
    const activityRaw = this.activityColInput?.value ?? '(none)';
    const activityWeight = this.activityWeightInput.value ?? 0.3;
    const activitySuffix = (activityRaw && activityRaw !== '(none)' && activityWeight > 0) ?
      ` · activity-weighted by ${activityRaw} (weight ${fmt(activityWeight)})` :
      '';
    this.summaryLabel.textContent =
      `Will score ${candidateCount} candidate${candidateCount === '1' ? '' : 's'} ` +
      `against ${refName} using ${presetKey} preset ` +
      `(Tc ${fmt(tcMin)}–${fmt(tcMax)}, atom-ratio ≤ ${fmt(mcs)}, CATS ≥ ${fmt(cats)})` +
      `${activitySuffix}.`;
  }

  /** Picks a short, human-friendly display name for the reference row.
   *  Tries common name-like columns ("name", "title", "label", "id",
   *  "compound") in that order; falls back to "row N" if nothing
   *  matches. Truncates long values so the summary stays one line. */
  private _referenceDisplayName(
    table: DG.DataFrame, refIdx: number, molCol: DG.Column,
  ): string {
    const candidates = ['name', 'Name', 'title', 'Title', 'label', 'Label',
      'id', 'ID', 'Id', 'compound', 'Compound'];
    for (const name of candidates) {
      const c = table.col(name);
      if (!c || c === molCol) continue;
      const v = c.get(refIdx);
      if (typeof v === 'string' && v.length > 0)
        return v.length > 40 ? v.substring(0, 37) + '...' : v;
    }
    return `row ${refIdx}`;
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
    this._saveMarksToCache(); // drops the cache entry — empty set is treated as "no marks"
    if (this._currentRefSmiles) this._reRender(this._currentRefSmiles);
    else this._updateSelectionCaption();
  }

  private _updateSelectionCaption() {
    const n = this.selectedAtoms.size;
    // Local-mode warning in the summary depends on whether marks exist,
    // so every mark-state change has to refresh it. Hooked here because
    // `_updateSelectionCaption` is the single chokepoint every paint /
    // erase / clear handler calls when marks change.
    this._updateSummary();
    if (n === 0) {
      // Modifier-key hint and "leave empty for any global hop" guidance
      // live permanently in `modifierHintLabel` below the preview now —
      // this caption only shows the dynamic "N atoms marked" status.
      this.selectionCaption.textContent = '';
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
      this.modifierHintLabel,
      this.selectionCaption,
      this.clearSelectionBtn,
    ], {style: {paddingRight: '12px'}});

    const right = ui.divV([
      this.summaryLabel,
      this.tableInput,
      this.colInputRoot,
      this.refRowInput,
      ui.h3('Threshold preset',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.presetInput,
      ui.divText('Sets Tc range + MCS atom-ratio cap. The numeric inputs ' +
        'below remain editable; the preset just seeds them. The category ' +
        'labels (Local / Easy / Middle / Hard) are informed by the ' +
        'qualitative hop taxonomy of Sun-Tawa-Wallqvist 2012 (DDT 17:310), ' +
        'which describes the regimes by NAME but does NOT prescribe numeric ' +
        'cutoffs — the thresholds below are Datagrok defaults, except ' +
        'Hard\'s ≤0.4 atom-ratio which matches Maeda 2024.',
      {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),
      ui.h3('Pre-filter shortlist (ECFP4 Tanimoto)',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.tanimotoMinInput,
      this.tanimotoMaxInput,
      this.rowCountPreviewLabel,
      ui.h3('Hop criterion (always applied)',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.mcsMaxInput,
      ui.divText('ratio_atom = atoms(MCS) / atoms(reference). A row passes ' +
        'this criterion iff ratio_atom ≤ max.',
      {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),

      ui.h3('Optional flag conditions',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.useCatsInFlagInput,
      this.catsSimInput,
      this.useTcInFlagInput,
      ui.divText('Each unchecked condition is computed and shown but does not ' +
        'affect the Scaffold Hop boolean flag. Uncheck both to keep only the ' +
        'atom-ratio criterion (≤ 0.4) driving the flag.',
      {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),

      this.showReasonInput,

      ui.h3('Activity-aware re-rank (optional)',
        {style: {marginTop: '12px', marginBottom: '4px'}}),
      this.activityColInputRoot,
      this.activityWeightInput,
      this.imputeActivityInput,
      this.mmpSupportFloorInput,
      this.predictKnownInput,
      ui.divText('Pick a numeric column (pIC50, pKi, etc.) to multiply the ' +
        'score by an activity-proximity factor. Leave as "(none)" to keep the ' +
        'structural-only ranking. The proximity itself is added as a separate ' +
        'column for inspection. Tick "Predict missing activity (MMP)" to also ' +
        'fill in masked / unknown activities via MMP rules. The two inputs ' +
        'below it are MMP knobs — leave at defaults unless you know why.',
      {style: {fontSize: '11px', color: 'var(--grey-5)', marginTop: '4px'}}),
    ], {style: {minWidth: '340px'}, classes: 'ui-form'});

    return ui.divH([left, right], {style: {alignItems: 'flex-start'}});
  }

  public getParams(): ScaffoldHoppingParams {
    const table = this.tableInput.value!;
    const molecules = this.colInput.value!;
    const referenceRowIdx = this._getCurrentRefIdx();
    if (referenceRowIdx == null)
      throw new Error('No reference row selected');
    if (this.tanimotoMinInput.value! > this.tanimotoMaxInput.value!)
      throw new Error('ECFP4 Tc min must be ≤ ECFP4 Tc max');
    const activityRaw = this.activityColInput?.value ?? '(none)';
    const activityColumnName = activityRaw === '(none)' ? '' : activityRaw;
    // R-group decomposition is the load-bearing primitive for Local
    // mode's exact-core replacement; the other presets use ErG. Read
    // the current preset choice from the dropdown rather than from the
    // module-level cache so any in-session preset switch is honored.
    const presetLabel = this.presetInput.value;
    const presetKey = (Object.keys(PRESETS) as PresetKey[])
      .find((k) => PRESET_LABELS[k] === presetLabel) ?? _lastPresetKey;
    const useRGroupReplacement = (presetKey === 'Local');
    // Local mode without marked atoms degenerates into "every near-
    // neighbour is a hop" — clicking through the orange-tinted warning
    // produces a useless 200+-row result table. Block here so the user
    // either marks atoms or switches preset. The error bubbles through
    // the dialog's `onOK` catch into grok.shell.error and the dialog
    // stays open so the user can correct.
    if (useRGroupReplacement && this.selectedAtoms.size === 0) {
      throw new Error('Local preset requires marked atoms. ' +
        'Ctrl-hover atoms on the reference molecule to mark the replaceable ' +
        'region, or switch to a different preset.');
    }
    return {
      table, molecules, referenceRowIdx,
      tanimotoMin: this.tanimotoMinInput.value!,
      tanimotoMax: this.tanimotoMaxInput.value!,
      mcsRatioMax: this.mcsMaxInput.value!,
      minCatsSim: this.catsSimInput.value!,
      replaceableAtoms: [...this.selectedAtoms].sort((a, b) => a - b),
      useTcInFlag: this.useTcInFlagInput.value ?? true,
      useCatsInFlag: this.useCatsInFlagInput.value ?? true,
      activityColumnName,
      activityWeight: this.activityWeightInput.value ?? 0.3,
      showReason: this.showReasonInput.value ?? false,
      useRGroupReplacement,
      imputeActivity: this.imputeActivityInput.value ?? false,
      mmpSupportFloor: this.mmpSupportFloorInput.value ?? 10,
      predictKnown: this.predictKnownInput.value ?? false,
    };
  }
}
