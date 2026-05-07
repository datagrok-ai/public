// Protein-Ligand interaction (ProLIF) context panel + batch handler — Docking variant.
//
// Mirrors `BiostructureViewer/src/utils/prolif-panel.ts` (intentional
// duplication until the shared surface moves to `@datagrok-libraries/bio`,
// see PR review). The Docking version differs in:
//   * `runPlBatchForDataset` pre-fetches the receptor once before the loop
//     (Docking poses don't carry a full receptor — ligand and receptor are
//     in separate files paired by REMARK metadata)
//   * the per-row guard is `isApplicableAutodock` (presence of a "binding
//     energy" REMARK), not `hasNonWaterHetatm`
//
// All other logic (column layout, conditional column inclusion, WeakMap
// renderer registry, error HTML, claim-cursor pool) is identical.

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getReceptorData, RECEPTOR_NAME_RE} from './utils';

// Kept in sync with `BiostructureViewer/scripts/protein_ligand_interactions.py:SKIP_RESNAMES`.
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

function detectNonWaterHetatmInstances(pdbText: string): string[] {
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

/** Source DataFrame context the batch handler needs. For Docking this is
 *  always populated with the same column for both `pdbCol` and `ligandCol`,
 *  signaling the handler to pre-fetch a single receptor and treat each
 *  row's pose as the ligand. */
export interface ProlifBatchCtx {
  df: DG.DataFrame;
  pdbCol: DG.Column<string>;
  ligandCol?: DG.Column<string>;
}

/** Builds the context-panel widget. Identical to BSV's makeProlifWidget. */
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
        // See BSV/src/utils/prolif-panel.ts for the rationale on iframe.onload
        // + fixed 600px height (the sandbox attribute blocks the parent from
        // measuring the iframe's contentDocument).
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
// Per-grid registry — see BSV/src/utils/prolif-panel.ts for the rationale.
const _plCellRendererWired = new WeakMap<DG.Grid, true>();

function _plErrorHtml(err: unknown): string {
  const msg = err instanceof Error ? err.message : String(err);
  const safe = msg.replace(/&/g, '&amp;').replace(/</g, '&lt;');
  return `<p style="color:#c0392b;padding:8px;font-family:sans-serif;margin:0">PL Error: ${safe}</p>`;
}

const PL_COL_NAMES_FIXED = {
  diag: 'PL Diagram',
  interactions: 'PL Interactions',
} as const;

// HBDonor / HBAcceptor are kept separate (chemistry direction matters).
// Other symmetric pairs collapse via INT_CODES on the Python side.
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

const ROW_PENDING = 0;
const ROW_SKIP = 1;
const ROW_OK = 2;
const ROW_ERROR = 3;

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

/** Per-row PL batch with Docking-specific receptor pre-fetch + AutoDock-pose
 *  precondition. See BSV/src/utils/prolif-panel.ts for the column-layout
 *  strategy. The receptor is resolved once from the first valid pose's
 *  REMARK header; getReceptorData makes 3 RPCs per call so we pre-screen
 *  poses with the cheap regex match first. */
export async function runPlBatchForDataset(ctx: ProlifBatchCtx): Promise<void> {
  const {df, pdbCol, ligandCol} = ctx;
  const rowCount = df.rowCount;
  const {diagName, resolve} = _resolvePlColumnNames(df);

  const diagCol = df.columns.addNewString(diagName) as DG.Column<string>;
  const interactionsCol = df.columns.addNewString(
    resolve(PL_COL_NAMES_FIXED.interactions)) as DG.Column<string>;
  const dropFixedCols = () => {
    df.columns.remove(diagCol.name);
    df.columns.remove(interactionsCol.name);
  };

  const counts: {[key: string]: Int32Array} = {};
  for (const c of PL_INT_COLS)
    counts[c.key] = new Int32Array(rowCount);
  const rowState = new Uint8Array(rowCount);

  // Docking-specific: pre-fetch the receptor once (shared across all poses).
  // Pre-screen poses with a regex match (cheap, no I/O) before the file
  // fetch — getReceptorData makes 3 RPCs per call.
  let receptor: string | null = null;
  if (ligandCol != null) {
    let firstPoseWithReceptor: string | null = null;
    for (let i = 0; i < rowCount; i++) {
      const pose = pdbCol.get(i);
      if (!pose || !isApplicableAutodock(pose)) continue;
      if (RECEPTOR_NAME_RE.test(pose)) {
        firstPoseWithReceptor = pose;
        break;
      }
    }
    if (firstPoseWithReceptor == null) {
      grok.shell.warning(
        'Could not resolve receptor for PL batch — no usable poses found.');
      dropFixedCols();
      return;
    }
    try {
      const data = await getReceptorData(firstPoseWithReceptor);
      receptor = typeof data.data === 'string'
        ? data.data
        : new TextDecoder().decode(data.data as Uint8Array);
    } catch (err) {
      grok.shell.warning(`Could not fetch receptor for PL batch: ${
        err instanceof Error ? err.message : String(err)}`);
      dropFixedCols();
      return;
    }
  }

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
    const pose = pdbCol.get(i);
    if (!pose || (ligandCol != null && !isApplicableAutodock(pose))) {
      diagCol.set(i, '');
      interactionsCol.set(i, '');
      rowState[i] = ROW_SKIP;
      return;
    }
    try {
      const result = await grok.functions.call(
        'BiostructureViewer:ProteinLigandInteractionDiagram', {
          protein: receptor ?? pose,
          ligand: ligandCol != null ? pose : '',
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
    dropFixedCols();
    grok.shell.info('PL diagram batch cancelled.');
    return;
  }

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
    }
  }

  _wireHtmlColumnRenderer(df, diagName);
}

/** Same heuristic as `Docking:isApplicableAutodock` — extracted here so
 *  the batch loop and panel can call it without going through the
 *  Datagrok function registry. */
function isApplicableAutodock(molecule: string): boolean {
  return molecule.includes('binding energy');
}
