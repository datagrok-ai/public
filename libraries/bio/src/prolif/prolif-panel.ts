// Protein-Ligand interaction (ProLIF) context panel + batch handler.
//
// Used by both `BiostructureViewer` and `Docking` via `@datagrok-libraries/bio/src/prolif/prolif-panel`.
// Lives in the bio library (rather than duplicated per-package) per PR review
// (https://github.com/datagrok-ai/public/pull/3783) — the only per-package
// differences are:
//   * the per-row "applicable" gate (`hasNonWaterHetatm` vs `isApplicableAutodock`)
//   * the optional pre-flight step (Docking pre-fetches a single receptor)
//   * how the row's `(protein, ligand)` pair maps onto the script's two inputs
// Those are passed in via the `RunPlBatchOptions` callbacks below; everything
// else (UI breakdown, capture container, cache, column layout) is shared.

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
export const PROLIF_SKIP_RESNAMES: ReadonlySet<string> = new Set([
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

/** Heuristic: PDB has both protein (ATOM records) and a non-water HETATM
 *  ligand. Used as the precondition guard for both the panel and the
 *  per-row batch loop, and exposed via a `@grok.decorators.func()` wrapper
 *  in `BiostructureViewer/src/package.ts` so the Datagrok panel system can
 *  use it as a condition. */
export function hasNonWaterHetatm(molecule: string): boolean {
  if (!molecule) return false;
  // Skip AutoDock poses — they're handled by the Docking package's own panel
  // and don't carry the receptor in the cell value.
  if (molecule.includes('binding energy')) return false;
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
// Per-row HTML cache (shared via window for cross-package consistency)
// ---------------------------------------------------------------------------

// Stored on `window` instead of at module level because the BSV and Docking
// packages still each import this from their own copy of `@datagrok-libraries/bio`
// when running unlinked — module-level Maps would then be DIFFERENT instances
// even though they look identical, causing the panel widget in one package
// to see an empty cache for data written by the other's batch handler.
// `df.name` is the inner key (stable across JS wrapper recreations, unlike
// the `DG.DataFrame` wrapper itself). The server-rendered PNG lives in the
// column value (cell renderer), not in this cache.
interface PlGlobalCache {
  html: Map<string, Map<number, string>>;
}
declare global {
  interface Window {
    __prolifCache?: PlGlobalCache;
  }
}
function _getCache(): PlGlobalCache {
  if (window.__prolifCache == null)
    window.__prolifCache = {html: new Map()};
  return window.__prolifCache;
}

/** Read the cached LigNetwork HTML for a (DataFrame, rowIdx) — null if no
 *  batch has produced one yet, or the cache was lost across packages. */
export function getPlHtmlForRow(df: DG.DataFrame, rowIdx: number): string | undefined {
  return _getCache().html.get(df.name)?.get(rowIdx);
}

function _setHtmlForRow(df: DG.DataFrame, rowIdx: number, html: string): void {
  const c = _getCache().html;
  let m = c.get(df.name);
  if (m == null) { m = new Map(); c.set(df.name, m); }
  m.set(rowIdx, html);
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

/** Semantic type stamped on the `PL Diagram` column so a dedicated context
 *  panel (registered in `BiostructureViewer/src/package.ts` via
 *  `@grok.decorators.panel`) fires when the user clicks any of its cells. */
export const PL_DIAGRAM_SEM_TYPE = 'PL.LigNetwork';

/** Given a `PL Diagram` column name, find the matching `PL Interactions`
 *  column. Modern runs always use canonical names (`PL Diagram` paired with
 *  `PL Interactions`), but the suffix-matching here is kept defensively so
 *  DataFrames carrying `(N)`-suffixed columns from older code keep working
 *  until the user re-runs the batch (which drops them and produces canonical
 *  names). Returns null if no matching column exists. */
export function interactionsColForDiagram(df: DG.DataFrame, diagramColName: string): DG.Column<string> | null {
  const m = /^PL Diagram( \(\d+\))?$/.exec(diagramColName);
  if (!m) return null;
  const suffix = m[1] ?? '';
  return (df.col(`PL Interactions${suffix}`) ?? null) as DG.Column<string> | null;
}

// Conditional per-type-count columns. Each entry maps a TS-side key to the
// Python script's column name plus the grid display name. HBDonor /
// HBAcceptor are kept separate because the chemistry distinction matters and
// a residue can be both (e.g. SER hydroxyl on different groups). CationPi/PiCation,
// XBDonor/XBAcceptor, MetalDonor/MetalAcceptor still collapse via INT_CODES on
// the Python side.
const PL_INT_COLS: ReadonlyArray<{key: string; pythonKey: string; displayName: string}> = [
  {key: 'nHBondDonor',    pythonKey: 'n_hbond_donor',    displayName: 'PL #H-bond donor'},
  {key: 'nHBondAcceptor', pythonKey: 'n_hbond_acceptor', displayName: 'PL #H-bond acceptor'},
  {key: 'nHydrophobic',   pythonKey: 'n_hydrophobic',    displayName: 'PL #Hydrophobic'},
  {key: 'nPiStacking',    pythonKey: 'n_pistacking',     displayName: 'PL #Pi-stacking'},
  {key: 'nCationic',      pythonKey: 'n_cationic',       displayName: 'PL #Cationic'},
  {key: 'nAnionic',       pythonKey: 'n_anionic',        displayName: 'PL #Anionic'},
  {key: 'nCationPi',      pythonKey: 'n_cationpi',       displayName: 'PL #Cation-Pi'},
  {key: 'nXBond',         pythonKey: 'n_xbond',          displayName: 'PL #X-bonds'},
  {key: 'nMetal',         pythonKey: 'n_metal',          displayName: 'PL #Metal'},
  {key: 'nVdw',           pythonKey: 'n_vdw',            displayName: 'PL #VdW'},
  {key: 'nTotal',         pythonKey: 'n_total',          displayName: 'PL #Total'},
];

// Matches any PL output column produced by this batch handler — canonical
// names (`PL Diagram`, `PL Interactions`, `PL #H-bond donor`, ...) and any
// `(N)` suffixed variants left behind by older versions that appended on
// re-run. `interactionsColForDiagram` still does suffix-aware lookups so
// stale state in user DataFrames keeps working.
const PL_COLUMN_PATTERN = /^PL (?:Diagram|Interactions|#[A-Za-z][\w -]*)( \(\d+\))?$/;

// Wipe every PL column the batch handler would have created in a prior run.
// Called at the start of `runPlBatch` so clicking "Compute for whole dataset"
// twice overwrites the columns in-place instead of producing a parallel
// `PL Diagram (2)` / `PL Interactions (2)` set.
function _dropPriorPlColumns(df: DG.DataFrame): void {
  const toDrop = df.columns.toList()
    .filter((c) => PL_COLUMN_PATTERN.test(c.name))
    .map((c) => c.name);
  for (const name of toDrop) df.columns.remove(name);
  // Also flush the cached LigNetwork HTML for this DataFrame — otherwise
  // the post-batch panel would still show stale diagrams for rows that
  // happen to be skipped on the re-run (e.g. user filtered the data).
  _getCache().html.delete(df.name);
}

// Set grid row height + column width for the diagram after the batch
// finishes. The actual cell rendering is handled by `PlLigNetworkPngRenderer`
// (registered in BSV/package.ts), which picks up the column via semType.
function _adjustGridForDiagrams(df: DG.DataFrame, colName: string): void {
  const tv = grok.shell.getTableView(df.name);
  if (tv == null) return;
  const gc = tv.grid.columns.byName(colName);
  if (gc != null) gc.width = 220;
  tv.grid.setOptions({rowHeight: 130});
  try { tv.grid.invalidate(); } catch (_e) { /* no-op */ }
}

// ---------------------------------------------------------------------------
// Offscreen vis.js canvas capture
// ---------------------------------------------------------------------------
//
// Capture-iframe hosting strategy: **30x30 container with overflow:hidden,
// full-size 600x520 iframe inside.**
//
// Why this and not the obvious approaches:
//   * `opacity: 0` / `opacity: 0.01` — Chrome's "visually imperceptible"
//     heuristic kicks in and throttles script execution inside the iframe.
//     vis.js never finishes drawing.
//   * `position: fixed; left: -99999px` — element is outside the viewport
//     bounds, Chrome's intersection-with-viewport heuristic flags it as
//     not visible, scripts throttle. vis.js never inserts the canvas.
//   * `visibility: hidden` / `display: none` — both block rendering
//     entirely; scripts may run but no canvas backing buffer exists.
//
// A 30x30 container + overflow:hidden works because the iframe's bounding
// rect is AT (0,0) with full 600x520 size, so Chrome's "in viewport" check
// passes. The parent container clips everything except a small visible
// region. Users see a small flickering thumbnail in the corner during a
// batch; auto-cleanup is via iframe.remove() in finish().
//
// No sandbox attribute on the iframe — we control the HTML (it's our own
// injected ProLIF output), so postMessage works and `iframe.contentDocument`
// is readable for the polling fallback.
//
// Capture concurrency = 3. Empirically, a SINGLE iframe in the small visible
// region gets throttled by Chrome's "visually imperceptible" heuristic
// (vis.js produces only a blank 10.3KB canvas). With multiple concurrent
// iframes Chrome stops throttling — each renders normally, and the 11KB-byte
// threshold below filters any blank captures that slip through.
const PL_CAPTURE_CONCURRENCY = 3;

let _captureContainer: HTMLDivElement | null = null;
function _getCaptureContainer(): HTMLDivElement {
  if (_captureContainer == null) {
    // Style lives in BSV's `viewer.css:.bsv-pl-capture-container`.
    _captureContainer = ui.div([], 'bsv-pl-capture-container');
    _captureContainer.setAttribute('data-prolif-capture', 'true');
    // `document.body.appendChild` is the one DOM operation without a clean
    // Datagrok wrapper: the capture container must be a `position:fixed`
    // overlay that survives every view switch (its lifetime is the whole
    // browser session, not any one TableView). `grok.shell.v.root` would
    // scope it to the current view and lose iframes mid-batch on view change.
    document.body.appendChild(_captureContainer);
  }
  return _captureContainer;
}

interface CaptureTask {
  df: DG.DataFrame;
  rowIdx: number;
  html: string;
  diagramCol: DG.Column<string>;
}
const _captureQueue: CaptureTask[] = [];
let _captureWorkersRunning = 0;

function _enqueueCapture(
  df: DG.DataFrame, rowIdx: number, html: string, diagramCol: DG.Column<string>,
): void {
  const existing = diagramCol.get(rowIdx);
  if (existing && existing.length > 100) return;
  _captureQueue.push({df, rowIdx, html, diagramCol});
  _kickCaptureWorkers();
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
      const {df, rowIdx, html, diagramCol} = task;
      const existing = diagramCol.get(rowIdx);
      if (existing && existing.length > 100) continue;
      const dataUrl = await _captureOne(html);
      if (dataUrl) {
        const base64 = dataUrl.replace(/^data:image\/png;base64,/, '');
        diagramCol.set(rowIdx, base64);
        const tv = grok.shell.getTableView(df.name);
        try { tv?.grid?.invalidate?.(); } catch (_e) { /* grid might be closing */ }
      }
    }
  } finally {
    _captureWorkersRunning--;
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
    const container = _getCaptureContainer();
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

/** Builds the ad-hoc context-panel widget shown when the user clicks a
 *  Molecule3D / PDB_ID / docking-pose cell. Just shows a brief ligand
 *  summary and the "Compute for whole dataset" button — the interactive
 *  LigNetwork diagram + breakdown live in the post-batch panel
 *  (`plDiagramInteractionsWidget` in BSV/package.ts), which fires on
 *  `PL Diagram` cell clicks.
 *
 *  `runBatch` overrides what the button does. The default runs `runPlBatch`
 *  with the BSV-style gate (`hasNonWaterHetatm`). The Docking package
 *  supplies its own runner that pre-fetches the receptor and gates rows
 *  with `isApplicableAutodock`. */
export function makeProlifWidget(params: {
  protein: string;
  ligand?: string;
  // Kept in signature for backward-compat with the panel decorator's
  // param shape; ignored now that we no longer compute on click.
  ligand_resname?: string;
}, batchCtx?: ProlifBatchCtx, runBatch?: () => void | Promise<void>): DG.Widget {
  const host = ui.div([], 'd4-empty-parent');

  const detectionSource = (params.ligand && params.ligand.trim()) || params.protein;
  const ligands = detectNonWaterHetatmInstances(detectionSource);

  // Style: `viewer.css:.bsv-pl-summary`.
  const summaryText = ligands.length === 0
    ? 'No non-water HETATM ligand found in this structure.'
    : ligands.length === 1
      ? `Ligand: ${ligands[0]}`
      : `${ligands.length} ligands found: ${ligands.join(', ')}`;
  host.append(ui.divText(summaryText, 'bsv-pl-summary'));

  const onClick = runBatch ?? (() => batchCtx != null
    ? runPlBatch({ctx: batchCtx, buildRowArgs: defaultBuildRowArgs(batchCtx)})
    : undefined);

  const btn = ui.button('Compute for whole dataset', onClick);
  btn.classList.add('bsv-pl-batch-btn');
  if (batchCtx == null && runBatch == null) {
    btn.disabled = true;
    ui.tooltip.bind(btn, 'No source DataFrame context — open this panel from a grid cell.');
  }
  host.append(btn);

  return new DG.Widget(host);
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

const PL_BATCH_CONCURRENCY = 4;

// Per-row state for the batch loop.
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

/** Per-row script calls via a 4-worker promise pool. We don't use the
 *  `meta:{vectorFunc:'true'}` Datagrok pattern (Chem's `getDescriptors`,
 *  `getMorganFingerprints`) because that's designed for column→column
 *  transforms with a single output value per row; ProLIF returns 13 values
 *  per row plus a large HTML blob, so a manual claim-cursor pool is the
 *  right fit.
 *
 *  Column layout:
 *    1. Add the `PL Diagram` and `PL Interactions` columns up front; the
 *       grid shows progress as workers complete.
 *    2. Accumulate per-type counts in Int32Arrays during the loop.
 *    3. After the loop, add only the count columns that have at least one
 *       OK row with value > 0 (skipping types with no observed instances —
 *       e.g. metal columns vanish for purely organic-ligand datasets).
 *
 *  Diagram columns are stamped with `PL_DIAGRAM_SEM_TYPE` so both
 *  `PlLigNetworkPngRenderer` (cell renderer) and `plDiagramInteractionsWidget`
 *  (post-batch panel) light up automatically. */
export async function runPlBatch(opts: RunPlBatchOptions): Promise<void> {
  const {ctx, buildRowArgs, prepare} = opts;
  const {df} = ctx;
  const rowCount = df.rowCount;

  // Drop any PL columns from a prior run so re-clicking the button overwrites
  // the canonical-name columns in place instead of producing a parallel
  // `PL Diagram (2)` / `PL Interactions (2)` set every click.
  _dropPriorPlColumns(df);

  // PL Diagram added BEFORE PL Interactions so it appears to the left in
  // the grid — the diagram is the primary visual artefact, the interaction
  // string is supporting detail.
  const diagramCol = df.columns.addNewString(PL_COL_NAMES_FIXED.diagram) as DG.Column<string>;
  diagramCol.semType = PL_DIAGRAM_SEM_TYPE;
  const interactionsCol = df.columns.addNewString(PL_COL_NAMES_FIXED.interactions) as DG.Column<string>;
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

  const counts: {[key: string]: Int32Array} = {};
  for (const c of PL_INT_COLS)
    counts[c.key] = new Int32Array(rowCount);
  const rowState = new Uint8Array(rowCount);

  const pi = DG.TaskBarProgressIndicator.create(
    `Computing PL interactions (${rowCount} rows)...`, {cancelable: true});
  let completed = 0;
  let cursor = 0;
  const claim = (): number => {
    if (pi.canceled) return -1;
    if (cursor >= rowCount) return -1;
    return cursor++;
  };

  const processRow = async (i: number): Promise<void> => {
    const args = buildRowArgs(i);
    if (args == null) {
      interactionsCol.set(i, '');
      rowState[i] = ROW_SKIP;
      return;
    }
    try {
      const result = await grok.functions.call(
        'BiostructureViewer:ProteinLigandInteractionDiagram', args,
      ) as DG.DataFrame;
      interactionsCol.set(i, result.col('interactions')!.get(0) as string);
      // Cache the LigNetwork HTML for the panel widget. Enqueue an offscreen
      // capture so vis.js renders, posts the canvas PNG back via postMessage,
      // and the captured base64 lands in `PL Diagram` for the cell renderer.
      const html = result.col('html')!.get(0) as string;
      _setHtmlForRow(df, i, html);
      _enqueueCapture(df, i, html, diagramCol);
      for (const c of PL_INT_COLS)
        counts[c.key][i] = result.col(c.pythonKey)!.get(0) as number;
      rowState[i] = ROW_OK;
    } catch (err) {
      interactionsCol.set(i, err instanceof Error ? `Error: ${err.message}` : String(err));
      rowState[i] = ROW_ERROR;
    }
  };

  try {
    await Promise.all(Array.from({length: PL_BATCH_CONCURRENCY}, async () => {
      for (;;) {
        const i = claim();
        if (i === -1) return;
        await processRow(i);
        completed++;
        pi.update(100 * completed / rowCount, `${completed} / ${rowCount}`);
      }
    }));
  } finally {
    pi.close();
  }

  if (pi.canceled) {
    dropFixedCols();
    grok.shell.info('PL interaction batch cancelled.');
    return;
  }

  // Add only the count columns that have at least one OK row with value > 0.
  // For each kept column, write counts for OK rows and DG.INT_NULL for ERROR
  // rows so per-column stats / sorts / filters distinguish "failed" from
  // "0 interactions". SKIP rows stay at the addNewInt default (DG.INT_NULL).
  for (const c of PL_INT_COLS) {
    const arr = counts[c.key];
    // Skip the whole column if no OK row has a positive count. `arr.some`
    // short-circuits like the manual loop+break but reads cleanly.
    if (!arr.some((v, i) => rowState[i] === ROW_OK && v > 0)) continue;
    const intCol = df.columns.addNewInt(c.displayName) as DG.Column<number>;
    for (let i = 0; i < rowCount; i++) {
      if (rowState[i] === ROW_OK) intCol.set(i, arr[i]);
      else if (rowState[i] === ROW_ERROR) intCol.set(i, DG.INT_NULL);
    }
  }

  _adjustGridForDiagrams(df, diagramCol.name);
  grok.shell.info(`PL interactions added for ${rowCount} rows.`);
}
