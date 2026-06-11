// Protein-Ligand interaction (ProLIF) context panel + batch handler.
//
// Used by:
//   * `pdbInteractionsWidget` panel — regular PDB Molecule3D cells
//   * `pdbIdInteractionsWidget` panel — PDB_ID strings, fetches from RCSB
//   * `PlDiagramObjectHandler` — `PL Diagram` cells, registered in `init()`
//     and keyed off the `.%prolif-source` column tag
//   * `dockingInteractionsWidget` panel — AutoDock-pose cells, via the
//     `./docking-pose-prolif.ts` wrapper that supplies the receptor
//     pre-fetch + `isAutoDockPose` gate as `RunPlBatchOptions` callbacks

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// ---------------------------------------------------------------------------
// HETATM detection
// ---------------------------------------------------------------------------

// Kept in sync with `BiostructureViewer/scripts/protein_ligand_interactions.py:SKIP_RESNAMES`.
// Common biological cofactors (HEM, NAD, FAD, ATP, etc.) are listed because
// they outvote a small inhibitor on the most-common-HETATM heuristic and
// would otherwise silently steer ProLIF onto the wrong ligand.
const PROLIF_SKIP_RESNAMES: ReadonlySet<string> = new Set([
  // waters
  'HOH', 'WAT', 'H2O', 'D2O', 'DOD',
  // ions / metals
  'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'MN', 'FE', 'CU', 'NI',
  // crystallographic buffers / additives
  'SO4', 'PO4', 'NO3', 'ACT', 'CO3',
  'GOL', 'EDO', 'PEG', 'PG4', 'DMS', 'TRS', 'IMD', 'BME',
  // common biological cofactors
  'HEM', 'HEC', 'HEB', 'HEA',
  'NAD', 'NAI', 'NAP', 'NDP', 'NAH', 'NAJ',
  'FAD', 'FMN', 'FDA',
  'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',
  'COA', 'ACO', 'COO',
  'SAM', 'SAH',
  'PLP', 'PMP',
  'BTN',
  'CLA', 'CHL',
  'B12', 'COB', 'BCA',
]);

/** Single source of truth for "looks like an AutoDock pose" — the
 *  `binding energy` REMARK line is what AutoDock writes into every pose
 *  output. BSV uses this to *reject* such inputs from `hasNonWaterHetatm`
 *  (BSV's `dockingInteractionsWidget` handles them instead). */
export function isAutoDockPose(molecule: string): boolean {
  return molecule.includes('binding energy');
}

/** Heuristic: PDB has both protein (ATOM records) and a non-water HETATM
 *  ligand. Used as the precondition guard for both the panel and the
 *  per-row batch loop, and exposed via a `@grok.decorators.func()` wrapper
 *  in `BiostructureViewer/src/package.ts` so the Datagrok panel system can
 *  use it as a condition. */
export function hasNonWaterHetatm(molecule: string): boolean {
  if (!molecule) return false;
  // AutoDock poses are handled by BSV's `dockingInteractionsWidget` panel
  // (they don't carry a full receptor in the cell value; the script would fail).
  if (isAutoDockPose(molecule)) return false;
  if (!/^ATOM/m.test(molecule)) return false;
  return /^HETATM.{11}(?!HOH|WAT|H2O)/m.test(molecule);
}

/** Returns each unique non-water HETATM ligand instance as a "RESNAME CHAIN RESID"
 *  identifier (e.g. "STI A 401"). When the same resname appears at multiple sites
 *  or in multiple chains, each instance gets its own entry so the user can pick.
 *  Common biological cofactors (HEM, NAD, FAD, ATP, ...) are filtered via
 *  `PROLIF_SKIP_RESNAMES` — see the comment on that constant for why. */
export function detectNonWaterHetatmInstances(pdbText: string): string[] {
  if (!pdbText) return [];
  const seen = new Map<string, number>();
  for (const line of pdbText.split('\n')) {
    if (line.startsWith('HETATM') && line.length >= 26) {
      const rn = line.slice(17, 20).trim();
      const chain = (line[21] || '').trim() || 'A';
      const resid = line.slice(22, 26).trim();
      if (rn && !PROLIF_SKIP_RESNAMES.has(rn)) {
        const key = `${rn} ${chain} ${resid}`;
        seen.set(key, (seen.get(key) || 0) + 1);
      }
    }
  }
  return Array.from(seen.entries()).sort((a, b) => b[1] - a[1]).map(([k]) => k);
}

// ---------------------------------------------------------------------------
// Per-row HTML cache
// ---------------------------------------------------------------------------

// The raw interactive LigNetwork HTML is stashed on the `PL Diagram`
// column's `temp` slot (Datagrok-canonical column-scoped session cache —
// same pattern as Chem's `col.temp[ChemTemps.SUBSTRUCT_PROVIDERS]` and
// Peptides' `col.temp[colorCacheKey]`). The captured PNG lives in the
// column value (read by PowerGrid's `rawPng` renderer); only the raw HTML
// for the post-batch context panel lives here.
//
// Lifetime = column lifetime: when `_dropPriorPlColumns` removes the
// column on a re-batch, the cache vanishes with it. Doesn't survive
// project save/reopen (the user would re-click the batch button), but
// neither does the offscreen-iframe machinery that produced the HTML.
export const PL_HTML_TEMP_KEY = '.prolif-html-cache';

/** Read the cached LigNetwork HTML for a row of the given `PL Diagram`
 *  column. Returns `undefined` if no batch has populated this column. */
export function getPlHtmlForRow(diagramCol: DG.Column, rowIdx: number): string | undefined {
  return (diagramCol.temp[PL_HTML_TEMP_KEY] as Map<number, string> | undefined)?.get(rowIdx);
}

function _getHtmlCache(diagramCol: DG.Column): Map<number, string> {
  let m = diagramCol.temp[PL_HTML_TEMP_KEY] as Map<number, string> | undefined;
  if (m == null) {
    m = new Map<number, string>();
    diagramCol.temp[PL_HTML_TEMP_KEY] = m;
  }
  return m;
}

// ---------------------------------------------------------------------------
// Per-interaction-type residue breakdown UI
// ---------------------------------------------------------------------------

// Interaction-type code → human-readable label, mirroring `INT_CODES` in
// the Python script. The display order is chemically grouped (donor/acceptor
// adjacent; positive/negative adjacent; etc.).
const INT_CODE_LABELS: Readonly<Record<string, string>> = {
  HD: 'H-bond donors',
  HA: 'H-bond acceptors',
  HY: 'Hydrophobic',
  PI: 'Pi-stacking',
  CAT: 'Cationic',
  AN: 'Anionic',
  CPI: 'Cation-Pi',
  XB: 'Halogen bonds',
  M: 'Metal',
  VDW: 'Van der Waals',
};
const INT_CODE_ORDER: readonly string[] = [
  'HD', 'HA', 'HY', 'PI', 'CAT', 'AN', 'CPI', 'XB', 'M', 'VDW',
];

// Sort residues by numeric residue id (so GLY3 < GLN21 < TYR123). Plain
// `.sort()` is alphabetic and produces confusing orderings (GLN21 before GLY3).
function _sortResidues(residues: Iterable<string>): string[] {
  const resIdRe = /(\d+)/;
  return [...residues].sort((a, b) => {
    const ma = resIdRe.exec(a);
    const mb = resIdRe.exec(b);
    const na = ma ? parseInt(ma[1], 10) : 0;
    const nb = mb ? parseInt(mb[1], 10) : 0;
    if (na !== nb) return na - nb;
    return a.localeCompare(b);
  });
}

/** Parse the comma-separated interactions string the Python script produces
 *  (`"ASP40_CAT, ASP40_HD, GLY11_VDW, ..."` — each entry is
 *  `<RESNAME+RESID>_<CODE>`) and return a context-panel UI that groups
 *  residues by interaction type. All styling lives in
 *  `viewer.css:.bsv-pl-interaction-breakdown-*`. */
export function renderInteractionBreakdown(interactionsStr: string): HTMLElement {
  const host = ui.div([], 'bsv-pl-interaction-breakdown');

  // Bucket residues by code. Each entry is `<residue>_<code>` — split at the
  // LAST underscore so residue names with internal digits/letters survive.
  const buckets: Record<string, Set<string>> = {};
  for (const raw of (interactionsStr || '').split(',')) {
    const entry = raw.trim();
    const i = entry.lastIndexOf('_');
    if (i < 0) continue;
    (buckets[entry.slice(i + 1)] ??= new Set()).add(entry.slice(0, i));
  }
  const totalCount = Object.values(buckets).reduce((s, set) => s + set.size, 0);

  host.append(ui.divH([
    ui.divText('Interactions', 'bsv-pl-interaction-breakdown-title'),
    ui.divText(`(${totalCount})`, 'bsv-pl-interaction-breakdown-count'),
  ], 'bsv-pl-interaction-breakdown-header'));

  if (totalCount === 0) {
    host.append(ui.divText('No interactions detected.', 'bsv-pl-interaction-breakdown-empty'));
    return host;
  }

  const addRow = (label: string, residues: Set<string>) => host.append(ui.divH([
    ui.divText(`${label} (${residues.size}):`, 'bsv-pl-interaction-breakdown-label'),
    ui.divText(_sortResidues(residues).join(', '), 'bsv-pl-interaction-breakdown-list'),
  ], 'bsv-pl-interaction-breakdown-row'));

  // Known codes in their canonical display order, then any unknown codes
  // (defensive — surface them under raw names rather than silently dropping).
  for (const code of INT_CODE_ORDER)
    if (buckets[code]?.size) addRow(INT_CODE_LABELS[code], buckets[code]);
  for (const code of Object.keys(buckets))
    if (!INT_CODE_ORDER.includes(code)) addRow(code, buckets[code]);
  return host;
}

// ---------------------------------------------------------------------------
// Column layout + semType
// ---------------------------------------------------------------------------

const PL_COL_NAMES_FIXED = {
  interactions: 'PL Interactions',
  diagram: 'PL Diagram',
} as const;

/** Semantic type stamped on the `PL Diagram` column. PowerGrid's built-in
 *  `RawPNGRenderer` (registered for `cellType: 'rawPng'`) picks columns up
 *  via this semType and draws the base64 PNG values directly — we no
 *  longer ship a custom cell renderer. The context panel is wired
 *  separately through `DG.ObjectHandler.register(new PlDiagramObjectHandler())`
 *  in BSV/package.ts's `init`, keyed off the `.%prolif-source` column tag. */
export const PL_DIAGRAM_SEM_TYPE = 'rawPng';

/** Tag set on every `PL Diagram` column by the batch handler. Its value is
 *  the name of the source ligand/PDB column the diagrams were computed
 *  from. Also used by `findPlDiagramColForSource` so the inline-diagram
 *  panel for a Molecule3D / docking-pose cell can hit the batch cache
 *  instead of re-running the script. */
export const PROLIF_SOURCE_TAG = '.%prolif-source';

/** Tag set on every `PL Diagram` column by the batch handler. Its value is
 *  the name of the sibling `PL Interactions` column. `interactionsColForDiagram`
 *  reads this tag rather than guessing the name by regex so a future
 *  rename (e.g. uniquification on column-name collision) doesn't break the
 *  pairing. */
export const PROLIF_INTERACTIONS_TAG = '.%prolif-interactions-col';

/** Look up the `PL Interactions` column paired with the given diagram column.
 *  Reads the `PROLIF_INTERACTIONS_TAG` tag set by `runPlBatch`. Returns null
 *  if the tag is missing (column wasn't produced by this PR's batch) or
 *  the referenced column has since been removed. */
export function interactionsColForDiagram(diagramCol: DG.Column): DG.Column<string> | null {
  const interactionsName = diagramCol.tags[PROLIF_INTERACTIONS_TAG];
  if (!interactionsName) return null;
  const df = diagramCol.dataFrame;
  if (df == null) return null;
  return (df.col(interactionsName) ?? null) as DG.Column<string> | null;
}

/** Find a `PL Diagram` column whose `.%prolif-source` tag points at the
 *  given source column. Returns null if no such column exists (no batch
 *  has been run for this source). Used by the inline-diagram path so
 *  clicking a Molecule3D / docking-pose cell can reuse the cached HTML
 *  from a prior batch instead of spending ~15-30s on a fresh script call. */
export function findPlDiagramColForSource(df: DG.DataFrame, sourceColName: string): DG.Column<string> | null {
  for (const col of df.columns.toList()) {
    if (col.semType === PL_DIAGRAM_SEM_TYPE && col.tags[PROLIF_SOURCE_TAG] === sourceColName)
      return col as DG.Column<string>;
  }
  return null;
}

// Conditional per-type-count columns. Each entry maps TS key → Python
// script column → grid display name. HBDonor/HBAcceptor stay separate
// because the chemistry distinction matters (a residue can be both);
// other symmetric pairs collapse on the Python side via INT_CODES. Short
// display codes keep grid columns narrow; full names appear in the
// breakdown UI where there's room.
const PL_INT_COLS: ReadonlyArray<{key: string; pythonKey: string; displayName: string}> = [
  {key: 'nHBondDonor', pythonKey: 'n_hbond_donor', displayName: 'PL HBD'},
  {key: 'nHBondAcceptor', pythonKey: 'n_hbond_acceptor', displayName: 'PL HBA'},
  {key: 'nHydrophobic', pythonKey: 'n_hydrophobic', displayName: 'PL Hyd'},
  {key: 'nPiStacking', pythonKey: 'n_pistacking', displayName: 'PL Pi'},
  {key: 'nCationic', pythonKey: 'n_cationic', displayName: 'PL Cat'},
  {key: 'nAnionic', pythonKey: 'n_anionic', displayName: 'PL An'},
  {key: 'nCationPi', pythonKey: 'n_cationpi', displayName: 'PL CatPi'},
  {key: 'nXBond', pythonKey: 'n_xbond', displayName: 'PL XB'},
  {key: 'nMetal', pythonKey: 'n_metal', displayName: 'PL Metal'},
  {key: 'nVdw', pythonKey: 'n_vdw', displayName: 'PL vdW'},
  {key: 'nTotal', pythonKey: 'n_total', displayName: 'PL Total'},
];

// Matches any PL output column we've ever produced — canonical names
// (`PL Diagram`, `PL Interactions`, `PL HBD` from `PL_INT_COLS`), legacy
// long-form names (`PL #H-bond donor`) from older versions, and any
// `(N)`-suffixed variants. `_dropPriorPlColumns` uses this to clean up
// stale state on a re-batch.
const _PL_COUNT_DISPLAY_NAMES = PL_INT_COLS.map((c) => c.displayName.replace(/^PL /, ''));
const PL_COLUMN_PATTERN = new RegExp(
  `^PL (?:Diagram|Interactions|${_PL_COUNT_DISPLAY_NAMES.join('|')}|#[A-Za-z][\\w -]*)( \\(\\d+\\))?$`,
);

// Wipe every PL column the batch handler would have created in a prior run.
// Called at the start of `runPlBatch` so clicking "Compute for whole dataset"
// twice overwrites the columns in-place instead of producing a parallel
// `PL Diagram (2)` / `PL Interactions (2)` set.
function _dropPriorPlColumns(df: DG.DataFrame): void {
  const toDrop = df.columns.toList()
    .filter((c) => PL_COLUMN_PATTERN.test(c.name))
    .map((c) => c.name);
  // Removing the column drops its `temp` slot too, so the prior batch's
  // cached LigNetwork HTML disappears with the `PL Diagram` column.
  for (const name of toDrop) df.columns.remove(name);
}

// Set grid row height + column width for the diagram after the batch
// finishes. The actual cell rendering is handled by PowerGrid's built-in
// `RawPNGRenderer` (semType `rawPng`), registered platform-wide via PowerGrid.
// Look up the TableView by exact DataFrame reference rather than by name —
// `grok.shell.getTableView(df.name)` returns the wrong view if multiple
// TableViews share a name (the platform allows duplicates).
function _adjustGridForDiagrams(df: DG.DataFrame, colNames: string[]): void {
  const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame?.id === df?.id);
  if (tv == null) return;
  for (const colName of colNames) {
    const gc = tv.grid.columns.byName(colName);
    if (gc != null) gc.width = 250;
  }
  // tv.grid.setOptions({rowHeight: 130});
  // try { tv.grid.invalidate(); } catch (_e) { /* no-op */ }
}

// ---------------------------------------------------------------------------
// Offscreen vis.js canvas capture
// ---------------------------------------------------------------------------
//
// Strategy: 30x30 container with `overflow:hidden`, full-size 600x520 iframe
// inside. Chrome throttles "visually imperceptible" iframes — `opacity:0`,
// off-viewport positioning, `visibility:hidden`, `display:none` all stop
// vis.js from drawing. A small but on-screen container slips past the
// heuristic; the full-size iframe is clipped to 30px but vis.js still draws
// at production resolution. Users see a brief thumbnail flicker during a
// batch. Concurrency = 2 keeps at least one iframe rendering at any
// moment during a busy batch (a SOLO iframe in the small region still
// gets throttled empirically). The 11KB threshold on the polled canvas
// rejects blank captures.
// No sandbox — we control the HTML, so postMessage + canvas access work.
const PL_CAPTURE_CONCURRENCY = 2;

// Lazy + ref-counted container. Created on first enqueue, removed once
// the queue drains and all workers exit so the bottom-right thumbnail
// flicker isn't visible during idle.
let _captureContainer: HTMLDivElement | null = null;
function _ensureCaptureContainer(): HTMLDivElement {
  if (_captureContainer == null) {
    _captureContainer = ui.div([], 'bsv-pl-capture-container');
    _captureContainer.setAttribute('data-prolif-capture', 'true');
    // `position:fixed` overlay appended to `document.body` so it isn't
    // scoped to the current view (mid-batch view switches would
    // otherwise detach the iframes).
    document.body.appendChild(_captureContainer);
  }
  return _captureContainer;
}
function _disposeCaptureContainerIfIdle(): void {
  if (_captureContainer != null && _captureQueue.length === 0 && _captureWorkersRunning === 0) {
    try { _captureContainer.remove(); } catch (_e) { /* */ }
    _captureContainer = null;
  }
}

interface CaptureTask {
  html: string;
  resolve: (base64: string | null) => void;
}
const _captureQueue: CaptureTask[] = [];
let _captureWorkersRunning = 0;

/** Queue an offscreen vis.js capture and return a Promise that resolves
 *  with the base64-encoded PNG (or null on timeout / malformed canvas). */
function _enqueueCapture(html: string): Promise<string | null> {
  return new Promise<string | null>((resolve) => {
    _captureQueue.push({html, resolve});
    _kickCaptureWorkers();
  });
}

function _kickCaptureWorkers(): void {
  while (_captureWorkersRunning < PL_CAPTURE_CONCURRENCY && _captureQueue.length > 0)
    void _runCaptureWorker();
}

async function _runCaptureWorker(): Promise<void> {
  _captureWorkersRunning++;
  try {
    while (_captureQueue.length > 0) {
      const task = _captureQueue.shift()!;
      const dataUrl = await _captureOne(task.html);
      task.resolve(dataUrl == null ? null : dataUrl.replace(/^data:image\/png;base64,/, ''));
    }
  } finally {
    _captureWorkersRunning--;
    _disposeCaptureContainerIfIdle();
  }
}

/** Recursively search `doc` and its nested iframes for the first canvas.
 *  ProLIF's `plot_lignetwork` HTML wraps the vis.js network in a child
 *  `<iframe srcdoc=...>` (so the LigNetwork CSS doesn't leak into the
 *  surrounding page). Our capture iframe is the OUTER frame, but the
 *  canvas vis.js creates lives one level deeper. srcdoc iframes inherit
 *  the parent origin, so `iframe.contentDocument` is readable in the
 *  same-origin chain. */
function _findCanvas(doc: Document | null): HTMLCanvasElement | null {
  if (doc == null) return null;
  const direct = Array.from(doc.getElementsByTagName('canvas'))
    .find((c) => c.width > 0 && c.height > 0);
  if (direct) return direct;
  for (const frame of Array.from(doc.getElementsByTagName('iframe'))) {
    try {
      const found = _findCanvas(frame.contentDocument);
      if (found) return found;
    } catch (_e) { /* cross-origin frame — skip */ }
  }
  return null;
}

function _captureOne(html: string, timeoutMs = 15000): Promise<string | null> {
  return new Promise((resolve) => {
    const container = _ensureCaptureContainer();
    const iframe = ui.element('iframe') as HTMLIFrameElement;
    iframe.classList.add('bsv-pl-capture-iframe');
    iframe.srcdoc = html;
    // NO sandbox — same-origin iframe so postMessage + canvas access work.
    let resolved = false;
    let timer: number | null = null;
    let pollHandle: number | null = null;
    let stableCount = 0;
    let lastSize = 0;
    const finish = (png: string | null) => {
      if (resolved) return;
      resolved = true;
      window.removeEventListener('message', onMsg);
      if (timer != null) window.clearTimeout(timer);
      if (pollHandle != null) window.clearInterval(pollHandle);
      try { iframe.remove(); } catch (_e) { /* */ }
      resolve(png);
    };
    // Fast path — injected JS postMessages once vis.js stabilises. Usually
    // doesn't fire because ProLIF wraps the network in a nested srcdoc
    // iframe whose `network` var is function-local; polling below is the
    // real workhorse.
    const onMsg = (ev: MessageEvent) => {
      if (resolved) return;
      const d = ev.data as {type?: string; pngDataUrl?: string} | undefined;
      if (!d || d.type !== 'prolif-ready') return;
      if (iframe.contentWindow != null && ev.source !== iframe.contentWindow) return;
      if (d.pngDataUrl && d.pngDataUrl.length > 100) finish(d.pngDataUrl);
      else finish(null);
    };
    window.addEventListener('message', onMsg);
    timer = window.setTimeout(() => finish(null), timeoutMs);
    container.appendChild(iframe);
    // Polling fallback (the actual primary path). Starts 500ms after mount
    // so the iframe has time to load vis-network from CDN before we start
    // hammering its canvas.
    window.setTimeout(() => {
      if (resolved) return;
      pollHandle = window.setInterval(() => {
        if (resolved) return;
        try {
          const doc = iframe.contentDocument;
          if (!doc) return;
          const canvas = _findCanvas(doc);
          if (canvas == null) return;
          if (canvas.width === 0 || canvas.height === 0) return;
          const dataUrl = canvas.toDataURL('image/png');
          if (!dataUrl || dataUrl.length < 100) return;
          // An empty 585x500 canvas serializes to ~10.3KB consistently;
          // real LigNetwork PNGs vary widely but always exceed 11KB once
          // any nodes/edges land. 11KB threshold rejects the blank-canvas
          // state while accepting even minimal real networks.
          if (dataUrl.length < 11000) {
            stableCount = 0;
            lastSize = 0;
            return;
          }
          // Wait for the dataURL to stay stable for 3 consecutive polls
          // (~1.5s of no change → vis.js physics settled).
          if (dataUrl.length === lastSize) {
            stableCount++;
            if (stableCount >= 3) finish(dataUrl);
          } else {
            lastSize = dataUrl.length;
            stableCount = 0;
          }
        } catch (_e) { /* transient iframe access errors are fine */ }
      }, 500);
    }, 500);
  });
}

// ---------------------------------------------------------------------------
// Ad-hoc context-panel widget
// ---------------------------------------------------------------------------

/** Source DataFrame context the batch handler needs. */
export interface ProlifBatchCtx {
  df: DG.DataFrame;
  pdbCol: DG.Column<string>;
  ligandCol?: DG.Column<string>;
}

/** Resolve a `ProlifBatchCtx` from a SemanticValue passed to a panel widget.
 *  Tries the cell's DataFrame + column first; falls back to the current
 *  TableView's first Molecule3D column when invoked outside a cell context
 *  (e.g. from the menu). The `cell.dart != null` guard matches the Bio
 *  package's pattern at `Bio/src/widgets/to-atomic-level-widget.ts:17` —
 *  without it, `cell.dataFrame` returns a wrapper around a null dart that
 *  would let bad data through.
 *
 *  Returns `undefined` when no usable Molecule3D column is found — the
 *  caller passes that straight through to `makeProlifWidget` which disables
 *  the batch button with a tooltip.
 *
 *  `asDockingPose=true` routes the resolved column as BOTH `pdbCol` (drives
 *  the per-row loop) AND `ligandCol` (signals the batch handler to
 *  pre-fetch the receptor and treat each row's value as the ligand). */
export function resolveBatchCtxFromSemValue(
  molecule: DG.SemanticValue, asDockingPose = false,
): ProlifBatchCtx | undefined {
  const cell = molecule.cell;
  const hasValidCell =
    cell != null && cell.dart != null && cell.dataFrame != null && cell.column != null;
  let df: DG.DataFrame | null = null;
  let pdbCol: DG.Column<string> | null = null;
  if (hasValidCell) {
    df = cell.dataFrame;
    pdbCol = cell.column as DG.Column<string>;
  } else
    return undefined;

  if (df == null || pdbCol == null) return undefined;
  return asDockingPose ? {df, pdbCol, ligandCol: pdbCol} : {df, pdbCol};
}

/** Builds the ad-hoc context-panel widget shown when the user clicks a
 *  Molecule3D / PDB_ID / docking-pose cell. Just shows a brief ligand
 *  summary and the "Compute for whole dataset" button — the interactive
 *  LigNetwork diagram + breakdown live in the post-batch panel
 *  (`PlDiagramObjectHandler` registered in BSV's `init`), which fires on
 *  `PL Diagram` cell clicks.
 *
 *  `runBatch` overrides what the button does. The default runs `runPlBatch`
 *  with the BSV-style gate (`hasNonWaterHetatm`). The Docking-flavour runner
 *  (`docking-pose-prolif.ts`) pre-fetches the receptor and gates rows
 *  with `isAutoDockPose`. */
export interface ProlifWidgetOptions {
  /** Override the batch button's click handler. The default runs `runPlBatch`
   *  with the BSV-style gate (`hasNonWaterHetatm`). The Docking-flavour runner
   *  in `docking-pose-prolif.ts` pre-fetches the receptor and gates rows with `isAutoDockPose`. */
  runBatch?: () => void | Promise<void>;
  /** When true, compute and mount the interactive LigNetwork iframe + the
   *  per-residue interaction breakdown inline on panel open (between the
   *  summary text and the batch button). Costs one Python script call per
   *  panel mount (~15-30s); a loader covers the wait. */
  showInlineDiagram?: boolean;
}

export function makeProlifWidget(params: {
  protein: string;
  ligand?: string;
  ligand_resname?: string;
}, batchCtx?: ProlifBatchCtx, opts?: ProlifWidgetOptions): DG.Widget {
  const host = ui.div([], 'd4-empty-parent');

  const detectionSource = (params.ligand && params.ligand.trim()) || params.protein;
  const ligands = detectNonWaterHetatmInstances(detectionSource);

  // Style: `viewer.css:.bsv-pl-summary`.
  const summaryText = ligands.length === 0 ?
    'No non-water HETATM ligand found in this structure.' :
    ligands.length === 1 ?
      `Ligand: ${ligands[0]}` :
      `${ligands.length} ligands found: ${ligands.join(', ')}`;
  host.append(ui.divText(summaryText, 'bsv-pl-summary'));

  // Optional inline LigNetwork: 0 ligands / no explicit ligand → nothing;
  // 1 ligand or caller-supplied → auto-compute; N>1 without selection →
  // show a "Ligand" picker and wait for the user.
  const hasExplicitLigand =
    (params.ligand != null && params.ligand.trim().length > 0) ||
    (params.ligand_resname != null && params.ligand_resname.trim().length > 0);
  if (opts?.showInlineDiagram && (ligands.length > 0 || hasExplicitLigand)) {
    const diagramHost = ui.div([], 'bsv-pl-inline-diagram');
    // Look for a cached diagram from a prior `runPlBatch`. When found, the
    // inline panel skips the 15-30s script call and renders from `column.temp`
    // (same cache the `PL Diagram` cell context panel reads). Only valid for
    // the single-ligand auto-detect path; the multi-ligand picker re-runs
    // the script for each user pick since the batch only computed one ligand.
    const cache = (batchCtx != null && (ligands.length <= 1 || hasExplicitLigand)) ?
      _findCachedDiagram(batchCtx.df, batchCtx.pdbCol.name) : undefined;
    if (ligands.length > 1 && !hasExplicitLigand) {
      // `nullable: true` + `value: null` forces an intentional pick rather
      // than guessing via the most-common-HETATM heuristic.
      const picker = ui.input.choice('Ligand', {
        value: null,
        items: ligands,
        nullable: true,
        onValueChanged: (v: string | null) => {
          if (v == null) return;
          ui.empty(diagramHost);
          void _mountInlineDiagram(diagramHost, {...params, ligand_resname: v});
        },
      });
      host.append(picker.root);
      host.append(diagramHost);
    } else {
      host.append(diagramHost);
      void _mountInlineDiagram(diagramHost, params, cache);
    }
  }

  const onClick = opts?.runBatch ?? (() => batchCtx != null ?
    runPlBatch({ctx: batchCtx, buildRowArgs: defaultBuildRowArgs(batchCtx)}) :
    undefined);

  const btn = ui.button('Compute for whole dataset', onClick);
  btn.classList.add('bsv-pl-batch-btn');
  if (batchCtx == null && opts?.runBatch == null) {
    btn.disabled = true;
    ui.tooltip.bind(btn, 'No source DataFrame context — open this panel from a grid cell.');
  }
  host.append(btn);

  return new DG.Widget(host);
}

/** Resolve the (diagramCol, rowIdx) pair if a prior batch has cached
 *  HTML for the current row of `sourceColName`. Uses `df.currentRowIdx`
 *  for the row (kept in sync by the platform with the clicked cell), and
 *  matches the diagram column via the `.%prolif-source` tag. Returns
 *  undefined if no batch has run or no HTML is cached for this row. */
function _findCachedDiagram(
  df: DG.DataFrame, sourceColName: string,
): {diagramCol: DG.Column<string>; rowIdx: number} | undefined {
  const rowIdx = df.currentRowIdx;
  if (rowIdx < 0) return undefined;
  const diagramCol = findPlDiagramColForSource(df, sourceColName);
  if (diagramCol == null) return undefined;
  if (getPlHtmlForRow(diagramCol, rowIdx) == null) return undefined;
  return {diagramCol, rowIdx};
}

/** Render the cached LigNetwork iframe + per-residue breakdown into `target`
 *  from the `(diagramCol, rowIdx)` pair. Pulls HTML from `column.temp`
 *  (populated by `runPlBatch`) and the interactions string from the
 *  tag-linked `PL Interactions` column. */
function _mountCachedDiagram(
  target: HTMLElement, diagramCol: DG.Column<string>, rowIdx: number,
): void {
  const html = getPlHtmlForRow(diagramCol, rowIdx)!;
  const iframe = ui.element('iframe') as HTMLIFrameElement;
  iframe.srcdoc = html;
  iframe.classList.add('bsv-pl-panel-iframe');
  iframe.onload = () => iframe.classList.add('bsv-pl-panel-iframe-loaded');
  target.append(iframe);
  const interactionsCol = interactionsColForDiagram(diagramCol);
  const interactionsStr = (interactionsCol?.get(rowIdx) as string | null) ?? '';
  target.append(renderInteractionBreakdown(interactionsStr));
}

/** Fetches the interactive LigNetwork + per-residue breakdown for a single
 *  (protein, ligand) pair via the Python script, then mounts them in the
 *  given target. Used by `makeProlifWidget`'s inline-diagram path. If `cache`
 *  is supplied, skips the script entirely and renders from `column.temp`
 *  (same source the `PL Diagram` cell context panel uses). The iframe is
 *  same-origin (no sandbox) so vis.js can read its own canvas; CSS lives
 *  in `viewer.css:.bsv-pl-panel-iframe`. */
async function _mountInlineDiagram(
  target: HTMLElement, params: {protein: string; ligand?: string; ligand_resname?: string},
  cache?: {diagramCol: DG.Column<string>; rowIdx: number},
): Promise<void> {
  if (cache != null) {
    _mountCachedDiagram(target, cache.diagramCol, cache.rowIdx);
    return;
  }
  const loader = ui.loader();
  target.append(loader);
  try {
    const result = await grok.functions.call(
      'BiostructureViewer:ProteinLigandInteractionDiagram', {
        protein: params.protein,
        ligand: params.ligand ?? '',
        ligand_resname: params.ligand_resname ?? '',
      },
    ) as DG.DataFrame;
    const html = result.col('html')!.get(0) as string;
    const interactionsStr = (result.col('interactions')?.get(0) as string | null) ?? '';
    const iframe = ui.element('iframe') as HTMLIFrameElement;
    iframe.srcdoc = html;
    iframe.classList.add('bsv-pl-panel-iframe');
    iframe.onload = () => iframe.classList.add('bsv-pl-panel-iframe-loaded');
    if (loader.isConnected) loader.remove();
    target.append(iframe);
    target.append(renderInteractionBreakdown(interactionsStr));
  } catch (err) {
    if (loader.isConnected) loader.remove();
    target.append(ui.divText(
      `Could not compute interactions: ${err instanceof Error ? err.message : String(err)}`));
  }
}

/** Default per-row argument builder. Skips rows that fail `hasNonWaterHetatm`. */
function defaultBuildRowArgs(ctx: ProlifBatchCtx): (rowIdx: number) => PlScriptArgs | null {
  return (rowIdx) => {
    const pdb = ctx.pdbCol.get(rowIdx);
    if (!pdb || !hasNonWaterHetatm(pdb)) return null;
    return {
      protein: pdb,
      ligand: ctx.ligandCol != null ? (ctx.ligandCol.get(rowIdx) ?? '') : '',
      ligand_resname: '',
    };
  };
}

// ---------------------------------------------------------------------------
// Per-row batch handler
// ---------------------------------------------------------------------------

// Script calls run sequentially. The Python script is single-threaded and
// the Datagrok server has its own worker scheduler; fanning out 4 concurrent
// calls from the client just queues them up server-side. Captures still
// pipeline against the next row's script call via `_enqueueCapture` (which
// is fire-and-forget) so wall-clock isn't worse than the previous design.

// Per-row state.
const ROW_SKIP = 1;
const ROW_OK = 2;
const ROW_ERROR = 3;

/** Arguments passed to `BiostructureViewer:ProteinLigandInteractionDiagram`. */
export interface PlScriptArgs {
  protein: string;
  ligand: string;
  ligand_resname: string;
}

/** Options for the per-row batch handler. `buildRowArgs` is the only
 *  required hook — it produces the script args for a row, or null to mark
 *  the row as skipped. `prepare` runs once before the loop and can return
 *  false to abort the whole run (e.g. Docking pre-fetches a single receptor
 *  and bails the batch if no usable pose is found). */
export interface RunPlBatchOptions {
  ctx: ProlifBatchCtx;
  buildRowArgs: (rowIdx: number) => PlScriptArgs | null;
  prepare?: () => Promise<boolean>;
}

/** Sequential per-row script calls (the Python script is single-threaded
 *  and the Datagrok server has its own worker scheduler, so client-side
 *  parallel fan-out doesn't help — see comment at the top of this section).
 *  Offscreen LigNetwork captures still pipeline against the next row's
 *  script call via the fire-and-forget `_enqueueCapture` queue.
 *
 *  We don't use `meta:{vectorFunc:'true'}` because the script returns 13
 *  values per row plus a large HTML blob plus a client-side capture step;
 *  a manual loop is the right fit.
 *
 *  Column layout:
 *    1. Add the `PL Diagram` and `PL Interactions` columns up front so the
 *       grid shows the structure immediately.
 *    2. Accumulate per-type counts in Int32Arrays during the loop.
 *    3. After the loop, bulk-fill all output columns via `col.init`
 *       (single notification per column instead of one per row).
 *
 *  Diagram columns are stamped with `PL_DIAGRAM_SEM_TYPE` (`rawPng`) so
 *  PowerGrid's built-in renderer draws the thumbnails and
 *  `PlDiagramObjectHandler` fires the context panel on cell click. */
export async function runPlBatch(opts: RunPlBatchOptions): Promise<void> {
  const {ctx, buildRowArgs, prepare} = opts;
  const {df} = ctx;
  const rowCount = df.rowCount;

  // Drop any PL columns from a prior run so re-clicking the button overwrites
  // the canonical-name columns in place instead of producing a parallel
  // `PL Diagram (2)` / `PL Interactions (2)` set every click.
  _dropPriorPlColumns(df);

  // PL Diagram added BEFORE PL Interactions so it appears to the left in
  // the grid. semType `rawPng` lights up PowerGrid's built-in PNG
  // renderer; the `.%prolif-source` tag links each cell to the source
  // ligand/PDB column (what `PlDiagramObjectHandler` and the inline-diagram
  // cache lookup key off); the `.%prolif-interactions-col` tag pairs the
  // diagram column with its sibling so `interactionsColForDiagram` doesn't
  // have to guess by name.
  const diagramCol = df.columns.addNewString(PL_COL_NAMES_FIXED.diagram) as DG.Column<string>;
  diagramCol.semType = PL_DIAGRAM_SEM_TYPE;
  diagramCol.tags[PROLIF_SOURCE_TAG] = ctx.pdbCol.name;
  const interactionsCol = df.columns.addNewString(PL_COL_NAMES_FIXED.interactions) as DG.Column<string>;
  diagramCol.tags[PROLIF_INTERACTIONS_TAG] = interactionsCol.name;
  // `<RESIDUE>_<CODE>` entries are joined with `, ` by the Python script.
  // Setting the multi-value separator lights up the platform's multi-value
  // filter mode so the user can filter for "engages GLY351" or "any HBDonor"
  // by ticking categories rather than typing substring queries.
  interactionsCol.tags[DG.TAGS.MULTI_VALUE_SEPARATOR] = ',';
  interactionsCol.tags[DG.TAGS.CELL_RENDERER] = 'Tags';
  interactionsCol.semType = 'Tags';
  const dropFixedCols = () => {
    df.columns.remove(interactionsCol.name);
    df.columns.remove(diagramCol.name);
  };

  // Optional pre-flight (e.g. Docking pre-fetches the receptor). Bail if it
  // signals abort.
  if (prepare != null) {
    let ok = false;
    try { ok = await prepare(); } catch (_e) { ok = false; }
    if (!ok) {
      dropFixedCols();
      return;
    }
  }

  // Per-row buffers. Column writes happen in a single bulk `col.init` pass
  // after the loop so the grid notifies on each output column exactly once.
  const interactionsValues = new Array<string>(rowCount).fill('');
  const counts: {[key: string]: Int32Array} = {};
  for (const c of PL_INT_COLS)
    counts[c.key] = new Int32Array(rowCount);
  const rowState = new Uint8Array(rowCount);
  const htmlByRow = _getHtmlCache(diagramCol);
  const pngByRow = new Map<number, string>();
  const capturePromises: Promise<void>[] = [];

  const pi = DG.TaskBarProgressIndicator.create(
    `Computing PL interactions (${rowCount} rows)...`, {cancelable: true});

  try {
    for (let i = 0; i < rowCount; i++) {
      if (pi.canceled) break;
      const args = buildRowArgs(i);
      if (args == null) {
        rowState[i] = ROW_SKIP;
        pi.update(100 * (i + 1) / rowCount, `${i + 1} / ${rowCount}`);
        continue;
      }
      try {
        const result = await grok.functions.call(
          'BiostructureViewer:ProteinLigandInteractionDiagram', args,
        ) as DG.DataFrame;
        interactionsValues[i] = result.col('interactions')!.get(0) as string;
        const html = result.col('html')!.get(0) as string;
        // Cache the interactive HTML on the diagram column's `temp` slot
        // for the post-batch context panel. Lifetime = column lifetime.
        htmlByRow.set(i, html);
        // Capture runs in the background and resolves after we've moved on
        // to the next row's script call — keeps the script→capture pipeline
        // saturated without parallel script calls.
        capturePromises.push(_enqueueCapture(html).then((base64) => {
          if (base64 != null) pngByRow.set(i, base64);
        }));
        for (const c of PL_INT_COLS)
          counts[c.key][i] = result.col(c.pythonKey)!.get(0) as number;
        rowState[i] = ROW_OK;
      } catch (err) {
        interactionsValues[i] = err instanceof Error ? `Error: ${err.message}` : String(err);
        rowState[i] = ROW_ERROR;
      }
      pi.update(100 * (i + 1) / rowCount, `${i + 1} / ${rowCount}`);
    }
  } finally {
    pi.close();
  }

  if (pi.canceled) {
    dropFixedCols();
    grok.shell.info('PL interaction batch cancelled.');
    return;
  }

  // Wait for the offscreen captures to drain — most have completed in
  // parallel with subsequent script calls, but the last few rows always
  // have outstanding captures at this point.
  if (capturePromises.length > 0) {
    const captureProgress = DG.TaskBarProgressIndicator.create(
      'Finalizing PL Diagram thumbnails...');
    try {
      await Promise.all(capturePromises);
    } finally {
      captureProgress.close();
    }
  }

  // Bulk-fill output columns. `col.init` writes the whole column without
  // firing per-row notifications — orders of magnitude fewer events than
  // a per-row `col.set` loop on a 1000-row dataset.
  interactionsCol.init((i) => interactionsValues[i]);
  diagramCol.init((i) => pngByRow.get(i) ?? null);
  // Add only the count columns that have at least one OK row with value > 0.
  // For each kept column, OK rows get their count, ERROR rows get
  // DG.INT_NULL so per-column stats / sorts / filters distinguish "failed"
  // from "0 interactions". SKIP rows also get DG.INT_NULL.
  for (const c of PL_INT_COLS) {
    const arr = counts[c.key];
    if (!arr.some((v, i) => rowState[i] === ROW_OK && v > 0)) continue;
    const intCol = df.columns.addNewInt(c.displayName) as DG.Column<number>;
    intCol.init((i) => rowState[i] === ROW_OK ? arr[i] : DG.INT_NULL);
  }

  _adjustGridForDiagrams(df, [diagramCol.name, interactionsCol.name]);
  grok.shell.info(`PL interactions added for ${rowCount} rows.`);
}
