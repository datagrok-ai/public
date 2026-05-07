// Protein-Ligand interaction (ProLIF) context panel + batch handler.
//
// Used from `package.ts` by both:
//   * the per-cell context-panel widget (`pdbInteractionsWidget`,
//     `pdbIdInteractionsWidget`)
//   * the per-row batch button on those widgets, which fans out the
//     `BiostructureViewer:ProteinLigandInteractionDiagram` script call
//     across rows of a DataFrame.
//
// Duplicated in `Docking/src/utils/prolif-panel.ts` until the shared
// surface lives in `@datagrok-libraries/bio` (we tried — see PR review).
// The Docking copy adds a receptor pre-fetch (the per-row pose value
// doesn't carry a full receptor; ligand+receptor have to be paired) but
// is otherwise byte-for-byte identical to this one.

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Kept in sync with `protein_ligand_interactions.py:SKIP_RESNAMES`. Common
// biological cofactors (HEM, NAD, FAD, ATP, etc.) are listed because they
// outvote a small inhibitor on the most-common-HETATM heuristic and would
// otherwise silently steer ProLIF onto the wrong ligand.
const PROLIF_SKIP_RESNAMES = new Set([
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
 *  in `package.ts` so the Datagrok panel system can use it as a condition. */
export function hasNonWaterHetatm(molecule: string): boolean {
  if (!molecule) return false;
  // Skip AutoDock poses — they're handled by the Docking package's own panel
  // and don't carry the receptor in the cell value.
  if (molecule.includes('binding energy')) return false;
  // Require both protein (ATOM records) and a non-water HETATM ligand —
  // otherwise the script has nothing meaningful to compute.
  if (!/^ATOM/m.test(molecule)) return false;
  return /^HETATM.{11}(?!HOH|WAT|H2O)/m.test(molecule);
}

/** Returns each unique non-water HETATM ligand instance as a "RESNAME CHAIN RESID"
 *  identifier (e.g. "STI A 401"). When the same resname appears at multiple sites
 *  or in multiple chains, each instance gets its own entry so the user can pick.
 *  Common biological cofactors (HEM, NAD, FAD, ATP, ...) are filtered via
 *  `PROLIF_SKIP_RESNAMES` — see the comment on that constant for why.
 *  Exported for testing; the panel widget itself uses it internally. */
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

/** Source DataFrame context the batch handler needs to write back into.
 *  `ligandCol` is optional: when present (Docking case), each row's pose
 *  becomes the ligand and a single shared receptor is the protein. When
 *  absent (BSV case), `pdbCol` carries protein + ligand in one PDB blob. */
export interface ProlifBatchCtx {
  df: DG.DataFrame;
  pdbCol: DG.Column<string>;
  ligandCol?: DG.Column<string>;
}

/** Builds the context-panel widget. Renders the diagram for the picked /
 *  inferred ligand, and adds a "Compute for whole dataset" batch button
 *  (disabled with tooltip when no batch context is wired up). */
export function makeProlifWidget(params: {
  protein: string;
  ligand?: string;
  ligand_resname?: string;
}, batchCtx?: ProlifBatchCtx): DG.Widget {
  const host = ui.div([], 'd4-empty-parent');
  const body = ui.div();
  host.append(body);

  const detectionSource = (params.ligand && params.ligand.trim()) || params.protein;
  const ligands = detectNonWaterHetatmInstances(detectionSource);

  const compute = (resname: string, target: HTMLElement) => {
    ui.empty(target);
    const loader = ui.loader();
    target.append(loader);
    (async () => {
      try {
        const result = await grok.functions.call(
          'BiostructureViewer:ProteinLigandInteractionDiagram',
          {protein: params.protein, ligand: params.ligand ?? '', ligand_resname: resname},
        ) as DG.DataFrame;
        const html = result.col('html')!.get(0) as string;
        const iframe = ui.element('iframe') as HTMLIFrameElement;
        iframe.srcdoc = html;
        iframe.classList.add('pl-panel-iframe');
        iframe.setAttribute('sandbox', 'allow-scripts');
        // Reveal once the iframe DOM has loaded — vis.js inside the iframe
        // finishes its async draw shortly after `load` fires. We use a
        // fixed 600px height because the sandbox attribute (without
        // `allow-same-origin`) blocks the parent from reading the iframe's
        // contentDocument to measure the rendered network size.
        const reveal = () => {
          if (loader.isConnected) loader.remove();
          iframe.style.opacity = '1';
        };
        iframe.onload = reveal;
        target.append(iframe);
      } catch (err) {
        ui.empty(target);
        target.append(ui.divText(
          `Could not compute interactions: ${err instanceof Error ? err.message : String(err)}`));
      }
    })();
  };

  if (ligands.length === 0 && !params.ligand_resname)
    body.append(ui.divText('No non-water HETATM ligand found in this structure.'));
  else if (params.ligand_resname || ligands.length === 1)
    compute(params.ligand_resname || ligands[0], body);
  else {
    const innerBody = ui.div();
    const picker = ui.input.choice('Ligand', {
      value: null, items: ligands, nullable: true,
      onValueChanged: (v: string | null) => { if (v) compute(v, innerBody); },
    });
    body.append(ui.divV([
      ui.divText(`${ligands.length} ligands found. Select one to compute interactions:`,
        {style: {marginBottom: '6px', color: 'var(--grey-5)'}}),
      picker.root, innerBody,
    ]));
  }

  // Always-rendered batch button. Disabled with tooltip when batchCtx missing.
  const btn = ui.button(
    'Compute for whole dataset',
    () => batchCtx != null ? runPlBatchForDataset(batchCtx) : undefined,
  );
  btn.classList.add('pl-batch-btn');
  if (batchCtx == null) {
    btn.setAttribute('disabled', '');
    btn.title = 'No source DataFrame context — open this panel from a grid cell.';
  }
  host.append(btn);

  return new DG.Widget(host);
}

const PL_BATCH_CONCURRENCY = 4;
// Per-grid registry of "PL cell renderer already wired here" so re-clicking
// the batch button doesn't stack listeners. WeakMap (not DataFrame.temp) so
// closing + reopening the table view (which destroys the old grid and creates
// a new one) gets a fresh registration on the new grid — without it, the
// reopened view would render PL Diagram cells as raw HTML strings.
const _plCellRendererWired = new WeakMap<DG.Grid, true>();

// The inline `style` is intentional — the cell renderer wraps this snippet
// in a sandboxed iframe via srcdoc, so its document doesn't inherit the
// parent's CSS and a class would do nothing. Escape `<` and `&` because the
// Python error message could contain markup.
function _plErrorHtml(err: unknown): string {
  const msg = err instanceof Error ? err.message : String(err);
  const safe = msg.replace(/&/g, '&amp;').replace(/</g, '&lt;');
  return `<p style="color:#c0392b;padding:8px;font-family:sans-serif;margin:0">PL Error: ${safe}</p>`;
}

// Always-present output columns: the diagram + the per-row interactions
// string. The per-type COUNT columns (PL #H-bond donor, etc.) are added
// conditionally — see PL_INT_COLS / runPlBatchForDataset below.
const PL_COL_NAMES_FIXED = {
  diag: 'PL Diagram',
  interactions: 'PL Interactions',
} as const;

// Conditional per-type-count columns. Each entry maps a TS-side key to:
//   * pythonKey: the column name in the Python script's result DataFrame
//   * displayName: the column name shown in the grid
// Order matters — it determines the column order in the grid output.
// HBDonor / HBAcceptor are kept separate because the chemistry distinction
// matters and a residue can be both (e.g. SER hydroxyl on different groups).
// CationPi/PiCation, XBDonor/XBAcceptor, MetalDonor/MetalAcceptor still
// collapse via INT_CODES on the Python side.
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

// Per-row state for the batch loop:
//   PENDING — row hasn't run yet (initial Uint8Array fill)
//   SKIP    — no ligand detected; row is empty
//   OK      — script returned a valid result; counts written to typed array
//   ERROR   — script threw; counts get DG.INT_NULL when columns are added
const ROW_PENDING = 0;
const ROW_SKIP = 1;
const ROW_OK = 2;
const ROW_ERROR = 3;

// Resolve a shared-suffix name for the diagram column ("PL Diagram",
// "PL Diagram (2)", ...) and pick the suffix every other column should
// share. If asymmetric prior state (user kept some, deleted others) would
// cause a name collision under the shared suffix, fall back to per-column
// getUnusedName — siblings then have mixed suffixes but won't collide.
function _resolvePlColumnNames(df: DG.DataFrame): {diagName: string; resolve: (base: string) => string} {
  const diagName = df.columns.getUnusedName(PL_COL_NAMES_FIXED.diag);
  const m = /^PL Diagram( \(\d+\))?$/.exec(diagName);
  const suffix = m && m[1] ? m[1] : '';
  const allBaseNames: string[] = [
    PL_COL_NAMES_FIXED.diag, PL_COL_NAMES_FIXED.interactions,
    ...PL_INT_COLS.map((c) => c.displayName),
  ];
  const sharedOk = allBaseNames.every(
    (base) => base === PL_COL_NAMES_FIXED.diag || df.col(`${base}${suffix}`) == null);
  const resolve = (base: string) =>
    sharedOk ? `${base}${suffix}` : df.columns.getUnusedName(base);
  return {diagName, resolve};
}

// Idempotent — sets `cellType='html'` + width on the named column, and
// installs ONE `onCellPrepare` subscriber per grid that renders every
// `PL Diagram*` column. Re-clicking the batch button on the same grid
// reuses the existing subscriber instead of stacking duplicates.
function _wireHtmlColumnRenderer(df: DG.DataFrame, colName: string) {
  const tv = grok.shell.getTableView(df.name);
  if (tv == null) {
    grok.shell.info(`Column "${colName}" added. Open a table view to see cell rendering.`);
    return;
  }
  const gc = tv.grid.columns.byName(colName);
  if (gc != null) {
    gc.cellType = 'html';
    gc.width = 600;
  }
  // Install the cell renderer subscription once per grid (not per column,
  // not per DataFrame), so re-clicks don't accumulate listeners and a
  // reopened view gets a fresh wire on the new grid instance. The handler
  // matches by name prefix — the `'PL Diagram'` prefix is part of this
  // column-set's contract: any column whose name starts with it will be
  // rendered as a sandboxed iframe (covers `PL Diagram`, `PL Diagram (2)`,
  // etc. — Datagrok's uniquify starts at 2).
  if (!_plCellRendererWired.has(tv.grid)) {
    tv.grid.onCellPrepare((cell) => {
      if (!cell.isTableCell || !cell.gridColumn.name.startsWith('PL Diagram')) return;
      const html = cell.cell.value as string | null;
      if (!html) return;
      const iframe = ui.element('iframe') as HTMLIFrameElement;
      iframe.srcdoc = html;
      iframe.className = 'pl-cell-iframe';
      iframe.setAttribute('sandbox', 'allow-scripts');
      cell.style.element = iframe;
    });
    _plCellRendererWired.set(tv.grid, true);
  }
  tv.grid.setOptions({'rowHeight': 550});
  grok.shell.info(`PL diagrams added to column "${colName}".`);
}

// Per-row script calls via a 4-worker promise pool. We don't use the
// `meta:{vectorFunc:'true'}` Datagrok pattern (Chem's `getDescriptors`,
// `getMorganFingerprints`) because that's designed for column→column
// transforms with a single output value per row; ProLIF returns 13 values
// per row plus a large HTML blob, so a manual claim-cursor pool is the
// right fit.
//
// Column layout strategy:
//   1. Add the "always-present" columns (diag + interactions) BEFORE the
//      loop so the cell renderer wires up immediately and rows visibly
//      populate as workers complete.
//   2. Accumulate per-type counts in Int32Arrays (one per type) during the
//      loop — no Datagrok column writes per type per row.
//   3. After the loop, ADD only the count columns that have at least one
//      OK row with value > 0 (skipping types with no observed instances —
//      e.g. metal columns vanish for purely organic-ligand datasets).
export async function runPlBatchForDataset(ctx: ProlifBatchCtx): Promise<void> {
  const {df, pdbCol, ligandCol} = ctx;
  const rowCount = df.rowCount;
  const {diagName, resolve} = _resolvePlColumnNames(df);

  const diagCol = df.columns.addNewString(diagName) as DG.Column<string>;
  const interactionsCol = df.columns.addNewString(
    resolve(PL_COL_NAMES_FIXED.interactions)) as DG.Column<string>;

  // Per-type accumulators. Keyed by PL_INT_COLS[i].key. Int32Array zero-fills.
  const counts: {[key: string]: Int32Array} = {};
  for (const c of PL_INT_COLS)
    counts[c.key] = new Int32Array(rowCount);
  // Per-row state (see ROW_* constants above).
  const rowState = new Uint8Array(rowCount);  // initialized to ROW_PENDING (0)

  const pi = DG.TaskBarProgressIndicator.create(
    `Computing PL diagrams (${rowCount} rows)...`, {cancelable: true});
  let completed = 0;
  let cursor = 0;
  const claim = (): number => {
    if (pi.canceled) return -1;
    if (cursor >= rowCount) return -1;
    return cursor++;
  };

  const processRow = async (i: number): Promise<void> => {
    const pdb = pdbCol.get(i);
    if (!pdb || !hasNonWaterHetatm(pdb)) {
      // No ligand to compute against — leave the string columns empty;
      // count columns will simply not get a value for this row (they'll
      // show empty, since addNewInt's default is DG.INT_NULL).
      diagCol.set(i, '');
      interactionsCol.set(i, '');
      rowState[i] = ROW_SKIP;
      return;
    }
    try {
      const result = await grok.functions.call(
        'BiostructureViewer:ProteinLigandInteractionDiagram', {
          protein: pdb,
          ligand: ligandCol != null ? (ligandCol.get(i) ?? '') : '',
          ligand_resname: '',
        },
      ) as DG.DataFrame;
      diagCol.set(i, result.col('html')!.get(0) as string);
      interactionsCol.set(i, result.col('interactions')!.get(0) as string);
      for (const c of PL_INT_COLS)
        counts[c.key][i] = result.col(c.pythonKey)!.get(0) as number;
      rowState[i] = ROW_OK;
    } catch (err) {
      diagCol.set(i, _plErrorHtml(err));
      interactionsCol.set(i, '');
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
    df.columns.remove(diagCol.name);
    df.columns.remove(interactionsCol.name);
    grok.shell.info('PL diagram batch cancelled.');
    return;
  }

  // Add only the count columns that have at least one OK row with value > 0.
  // For each kept column, write counts for OK rows and DG.INT_NULL for ERROR
  // rows so per-column stats / sorts / filters distinguish "failed" from
  // "0 interactions". SKIP and PENDING rows stay at the addNewInt default
  // (DG.INT_NULL, displayed as empty).
  for (const c of PL_INT_COLS) {
    const arr = counts[c.key];
    let hasAny = false;
    for (let i = 0; i < rowCount; i++) {
      if (rowState[i] === ROW_OK && arr[i] > 0) { hasAny = true; break; }
    }
    if (!hasAny) continue;
    const intCol = df.columns.addNewInt(resolve(c.displayName)) as DG.Column<number>;
    for (let i = 0; i < rowCount; i++) {
      if (rowState[i] === ROW_OK) intCol.set(i, arr[i]);
      else if (rowState[i] === ROW_ERROR) intCol.set(i, DG.INT_NULL);
      // ROW_SKIP / ROW_PENDING: leave at addNewInt default (DG.INT_NULL)
    }
  }

  _wireHtmlColumnRenderer(df, diagName);
}
