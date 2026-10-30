/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Platform} from '@datagrok-libraries/u2/dg';
import type {HitDesignApp} from '../hit-design-app';
import {HitAppBase} from '../hit-app-base';
import {HitDesignMergeConfig, HitDesignMergeMode, HitDesignTemplate, TemplateCompute} from '../types';
import {HitDesignMolColName, TileCategoriesColName, ViDColName, ViDSemType} from '../consts';
import {calculateCellValues} from './calculate-single-cell';
import {registerMolsBatch} from './molreg';
import {_package} from '../../package';

const MOL_NAME_HINTS = ['molecule', 'smiles', 'mol', 'structure', 'canonical_smiles', 'canonicalsmiles'];
const VID_NAME_HINTS = ['v-id', 'vid', 'v_id'];

// What a merge or a dropped file can be: csv, and every extension the platform has a reader for
// (HitAppBase.molFileExtReaders keeps only those).
export const MERGE_ALLOWED_EXTENSIONS: string[] = ['csv', ...HitAppBase.molFileExtReaders.map((r) => r.ext)];

export function getFileInfoExtension(fi: DG.FileInfo): string {
  const ext = fi.extension?.toLowerCase();
  if (ext) return ext;
  const name = fi.name?.toLowerCase() ?? '';
  const dot = name.lastIndexOf('.');
  return dot >= 0 ? name.slice(dot + 1) : '';
}

export async function readFileInfoToDataFrame(fi: DG.FileInfo): Promise<DG.DataFrame | null> {
  const ext = getFileInfoExtension(fi);
  return ext === 'csv' ? DG.DataFrame.fromCsv(await fi.readAsString()) : readMoleculeFile(ext, await fi.readAsBytes());
}

/** A file dropped from the user's computer. */
export async function readFileToDataFrame(file: File): Promise<DG.DataFrame | null> {
  const ext = file.name.split('.').pop()?.toLowerCase() ?? '';
  return ext === 'csv' ? DG.DataFrame.fromCsv(await file.text()) :
    readMoleculeFile(ext, new Uint8Array(await file.arrayBuffer()));
}

async function readMoleculeFile(ext: string, bytes: Uint8Array | null): Promise<DG.DataFrame | null> {
  const reader = HitAppBase.molFileExtReaders.find((r) => r.ext === ext);
  if (!reader || !bytes) return null;
  const bytesArg = reader.handlerFunc.inputs[0].name;
  const result: any = await reader.handlerFunc.apply({[bytesArg]: bytes});
  return Array.isArray(result) ? result[0] : result;
}

export function guessMoleculeColumn(cols: DG.Column[]): DG.Column | null {
  const semCols = cols.filter((c) => c.semType === DG.SEMTYPE.MOLECULE);
  if (semCols.length > 0) {
    const named = semCols.find((c) => MOL_NAME_HINTS.includes(c.name.toLowerCase()));
    return named ?? semCols[0];
  }
  const named = cols.find((c) => MOL_NAME_HINTS.includes(c.name.toLowerCase()));
  if (named) return named;
  return cols.find((c) => c.type === DG.COLUMN_TYPE.STRING) ?? null;
}

export function guessVidColumn(cols: DG.Column[]): DG.Column | null {
  const named = cols.find((c) => VID_NAME_HINTS.includes(c.name.toLowerCase()));
  if (named) return named;
  const semCol = cols.find((c) => c.semType === ViDSemType);
  if (semCol) return semCol;
  return null;
}

/** Marks a column the campaign view added empty for a layout's viewer bound to it: saved only once it holds values. */
export const PlaceholderTag = '.hit-design-placeholder';

/** A layout's column nothing has filled yet: a merge gives it the incoming column's type, a save leaves it out. */
export function isEmptyPlaceholder(col: DG.Column): boolean {
  return col.getTag(PlaceholderTag) === 'true' && col.stats.missingValueCount === col.length;
}

/** The campaign a merge goes into: the one the app has open, or the one an open is reading. */
export interface MergeTarget {
  df: DG.DataFrame;
  molColName: string;
  campaignId?: string;
  stages: string[];
  compute?: TemplateCompute;
}

function openTarget<T extends HitDesignTemplate>(app: HitDesignApp<T>): MergeTarget {
  if (!app.dataFrame) throw new Error('Campaign has no data');
  return {df: app.dataFrame, molColName: app.molColName, campaignId: app.campaignId, stages: app.stages,
    compute: app.campaign?.template?.compute ?? app.template?.compute};
}

export function isLocalUploadFileInfo(fi: DG.FileInfo): boolean {
  const fullPath = (fi.fullPath ?? '').toLowerCase();
  const name = (fi.name ?? '').toLowerCase();
  return !fullPath || fullPath === name;
}

type SourceMapping = {col: DG.Column; target: string; isComputed: boolean};

function getSourceMappings(
  incomingDf: DG.DataFrame,
  computedDf: DG.DataFrame | null,
  campMolColName: string,
  cfg: HitDesignMergeConfig,
): SourceMapping[] {
  const sources: SourceMapping[] = [];
  for (const col of Platform.listedColumns(incomingDf)) {
    // In SMILES mode, never carry over a V-iD-like column from the incoming table:
    // campaign V-iDs are authoritative.
    if (cfg.mode === 'smiles' && col.name === ViDColName) continue;

    let target = col.name;
    if (cfg.molColName && col.name === cfg.molColName) target = campMolColName;
    if (cfg.mode === 'vid' && cfg.vidColName && col.name === cfg.vidColName) target = ViDColName;

    sources.push({col, target, isComputed: false});
  }
  if (computedDf) {
    for (const col of Platform.listedColumns(computedDf)) {
      if (col.name === HitDesignMolColName || col.name === campMolColName) continue;
      if (sources.some((s) => s.target === col.name)) continue;
      sources.push({col, target: col.name, isComputed: true});
    }
  }
  return sources;
}

function setCorrectedValue(toCol: DG.Column, fromCol: DG.Column, toIdx: number, fromIdx: number): void {
  try {
    if (toCol.type === fromCol.type) {
      toCol.set(toIdx, fromCol.get(fromIdx), false);
      return;
    }
    if (toCol.type === DG.COLUMN_TYPE.STRING) {
      toCol.set(toIdx, fromCol.getString(fromIdx), false);
      return;
    }
    const v = fromCol.get(fromIdx);
    if (toCol.type === DG.COLUMN_TYPE.INT) {
      if (typeof v === 'number') toCol.set(toIdx, Math.round(v), false);
      else if (typeof v === 'string' && v.length > 0 && !isNaN(Number(v))) toCol.set(toIdx, Math.round(Number(v)), false);
      return;
    }
    if (toCol.type === DG.COLUMN_TYPE.FLOAT) {
      if (typeof v === 'number') toCol.set(toIdx, v, false);
      else if (typeof v === 'string' && v.length > 0 && !isNaN(Number(v))) toCol.set(toIdx, Number(v), false);
      return;
    }
    toCol.set(toIdx, v, false);
  } catch (e) {
    _package.logger.error(`Merge: failed to set ${toCol.name}[${toIdx}] from ${fromCol.name}[${fromIdx}]: ${e}`);
  }
}

function isCellEmpty(col: DG.Column, idx: number): boolean {
  if (col.isNone(idx)) return true;
  if (col.type === DG.COLUMN_TYPE.STRING) {
    const v = col.get(idx);
    return v === '' || v == null;
  }
  return false;
}

/**
 * Computes the campaign-side join keys, returns the per-row key array and a key → row map.
 * In SMILES mode the key is the canonical SMILES; in VID mode it's the V-iD string.
 */
function buildCampaignKeys(
  campDf: DG.DataFrame,
  campMolColName: string,
  mode: HitDesignMergeMode,
): {keys: (string | null)[]; map: Map<string, number>} {
  const keys: (string | null)[] = new Array(campDf.rowCount).fill(null);
  const map = new Map<string, number>();
  if (mode === 'smiles') {
    const col = campDf.col(campMolColName);
    if (!col) throw new Error(`Campaign has no molecule column "${campMolColName}"`);
    for (let i = 0; i < campDf.rowCount; i++) {
      const v = col.get(i);
      if (!v) continue;
      const canon = _package.convertToSmiles(v);
      if (canon && canon !== 'MALFORMED_INPUT_VALUE') {
        keys[i] = canon;
        if (!map.has(canon)) map.set(canon, i);
      }
    }
  } else {
    const col = campDf.col(ViDColName);
    if (!col) throw new Error(`Campaign has no "${ViDColName}" column`);
    for (let i = 0; i < campDf.rowCount; i++) {
      const v = col.get(i);
      if (!v) continue;
      keys[i] = v;
      if (!map.has(v)) map.set(v, i);
    }
  }
  return {keys, map};
}

function buildIncomingKeys(incomingDf: DG.DataFrame, cfg: HitDesignMergeConfig): (string | null)[] {
  const keys: (string | null)[] = new Array(incomingDf.rowCount).fill(null);
  if (cfg.mode === 'smiles') {
    const col = incomingDf.col(cfg.molColName!);
    if (!col) throw new Error(`Incoming table has no column "${cfg.molColName}"`);
    for (let i = 0; i < incomingDf.rowCount; i++) {
      const v = col.get(i);
      if (!v) continue;
      const canon = _package.convertToSmiles(v);
      if (canon && canon !== 'MALFORMED_INPUT_VALUE') keys[i] = canon;
    }
  } else {
    const col = incomingDf.col(cfg.vidColName!);
    if (!col) throw new Error(`Incoming table has no column "${cfg.vidColName}"`);
    for (let i = 0; i < incomingDf.rowCount; i++) {
      const v = col.get(i);
      if (v) keys[i] = String(v);
    }
  }
  return keys;
}

/**
 * Merges `incomingDf` into the current campaign of `app` (or `target`) based on the supplied config.
 * Supports both SMILES-based and V-iD-based joins, optional addition of new rows,
 * and optional execution of the campaign compute pipeline on the incoming molecules.
 */
export async function mergeIntoCampaign<T extends HitDesignTemplate>(
  app: HitDesignApp<T>,
  incomingDf: DG.DataFrame,
  cfg: HitDesignMergeConfig,
  target: MergeTarget = openTarget(app),
): Promise<void> {
  const campDf = target.df;
  const campMolColName = target.molColName;

  // Effective flags: silently strip combinations that the dialog should already have prevented.
  const canAddNewRows = cfg.addNewRows && !!cfg.molColName;
  const willRunCompute = cfg.runComputeOnNewRows && canAddNewRows;

  const {map: campKeyMap} = buildCampaignKeys(campDf, campMolColName, cfg.mode);
  const incKeys = buildIncomingKeys(incomingDf, cfg);

  // Decide which incoming rows match the campaign and which are net-new BEFORE running
  // compute, so the compute pipeline only runs on the rows that will actually be added.
  const matched: {campRow: number; incRow: number}[] = [];
  const unmatched: number[] = [];
  for (let i = 0; i < incomingDf.rowCount; i++) {
    const k = incKeys[i];
    if (!k) continue;
    const campRow = campKeyMap.get(k);
    if (campRow !== undefined) matched.push({campRow, incRow: i});
    else if (canAddNewRows) unmatched.push(i);
  }

  // Compute is only run for the unmatched (newly added) molecules, not for matched rows.
  // The resulting DataFrame is row-aligned with `unmatched` (index `u`), not with the
  // original incoming rows.
  let computedDf: DG.DataFrame | null = null;
  if (willRunCompute && cfg.molColName && unmatched.length > 0) {
    const compute = target.compute;
    const incMolCol = incomingDf.col(cfg.molColName);
    if (compute && incMolCol) {
      const newMols = unmatched.map((i) => incMolCol.get(i));
      const descriptors = compute.descriptors.enabled ? compute.descriptors.args ?? [] : [];
      try {
        computedDf = await calculateCellValues(newMols, descriptors, compute.functions, compute.scripts ?? [], compute.queries ?? []);
      } catch (e) {
        _package.logger.error(e);
        grok.shell.warning('Merge: compute pipeline failed; merging without computed columns.');
      }
    }
  }

  const sources = getSourceMappings(incomingDf, computedDf, campMolColName, cfg);

  // Add columns missing in the campaign; a layout's empty placeholder takes the incoming column's type
  const newlyAdded = new Set<string>();
  for (const s of sources) {
    const existing = campDf.col(s.target);
    if (existing !== null && !isEmptyPlaceholder(existing))
      continue;
    if (existing === null)
      campDf.columns.addNew(s.target, s.col.type);
    else if (existing.type !== s.col.type)
      campDf.changeColumnType(s.target, s.col.type as DG.COLUMN_TYPE);
    if (s.col.semType) campDf.col(s.target)!.semType = s.col.semType;
    newlyAdded.add(s.target);
  }

  const newRowStartIdx = campDf.rowCount;
  app.joining++;
  try {
    // Fill matched rows. Computed columns are skipped here — they were only computed for
    // the newly added subset and do not have values aligned with matched incoming rows.
    for (const {campRow, incRow} of matched) {
      for (const s of sources) {
        if (s.isComputed) continue;
        // Don't overwrite the join-key column on matched rows: it already matches by definition.
        if (cfg.mode === 'smiles' && s.target === campMolColName) continue;
        if (cfg.mode === 'vid' && s.target === ViDColName) continue;

        const targetCol = campDf.col(s.target);
        if (!targetCol) continue;
        const fillEmpty = isCellEmpty(targetCol, campRow);
        if (newlyAdded.has(s.target) || fillEmpty || cfg.clashStrategy === 'overwrite')
          setCorrectedValue(targetCol, s.col, campRow, incRow);
      }
    }

    // Add unmatched rows.
    const hasStages = target.stages.length > 0 && campDf.columns.contains(TileCategoriesColName);
    for (let u = 0; u < unmatched.length; u++) {
      const incRow = unmatched[u];
      campDf.rows.addNew(null, true);
      const newCampRow = campDf.rowCount - 1;
      for (const s of sources) {
        // Computed columns are aligned with the unmatched array (index `u`); incoming
        // columns stay aligned with the original incoming row index.
        const fromIdx = s.isComputed ? u : incRow;
        setCorrectedValue(campDf.col(s.target)!, s.col, newCampRow, fromIdx);
      }
      if (hasStages) {
        const stageCol = campDf.col(TileCategoriesColName)!;
        if (isCellEmpty(stageCol, newCampRow))
          stageCol.set(newCampRow, target.stages[0], false);
      }
    }
  } finally {
    app.joining--;
  }

  // Re-register new molecules so they get fresh, dictionary-backed V-iDs.
  // This deliberately overwrites any V-iD value carried over from the incoming table
  // (incoming V-iDs are not authoritative for the campaign).
  const incMolCol = cfg.molColName ? incomingDf.col(cfg.molColName) : null;
  const campVidCol = campDf.col(ViDColName);
  const newMols = incMolCol ? unmatched.map((i) => incMolCol.get(i)).filter((m) => !!m && m !== '') : [];
  if (incMolCol && campVidCol && target.campaignId && newMols.length > 0) {
    try {
      const vidMap = await registerMolsBatch(newMols, target.campaignId, app.appName);
      for (let j = 0; j < unmatched.length; j++) {
        const mol = incMolCol.get(unmatched[j]);
        const vid = mol ? vidMap.get(mol) : undefined;
        if (vid) campVidCol.set(newRowStartIdx + j, vid, false);
      }
    }
    catch (e) {
      _package.logger.error(e);
    }
  }

  campDf.fireValuesChanged();
}

/**
 * A save over a campaign another session saved since this one read it: the server's table is merged into `ours` by
 * canonical SMILES. Its rows we lack are added and its new columns copied into the rows we share, as before; in a
 * shared row a cell we left as we read it ({@link base}) takes the server's value, so the other session's edits are
 * kept, and a cell we edited keeps ours; a column we did not hold when we read the campaign we read as empty. Rows
 * without a valid molecule are never matched. Of a molecule more than one row holds on either side, which row is which
 * is a guess: its cells stay ours, and a new column takes the server's later row's value.
 */
export function mergeServerTable<T extends HitDesignTemplate>(app: HitDesignApp<T>, ours: DG.DataFrame,
  molColName: string, server: DG.DataFrame, base: DG.DataFrame | undefined): void {
  if (!server.col(molColName) || !ours.col(molColName))
    throw new Error('Molecule column not found');
  const key = (v: string | null): string | null => {
    const smiles = v ? _package.convertToSmiles(v) : null;
    return smiles && smiles !== 'MALFORMED_INPUT_VALUE' ? smiles : null;
  };
  /** Each row's key, the (last) row of each key, and the keys more than one row holds. */
  const rowsOf = (df: DG.DataFrame | undefined) => {
    const mols = df?.col(molColName);
    const keys = mols ? Array.from({length: df!.rowCount}, (_, i) => key(mols.get(i))) : [];
    const rows = new Map<string, number>();
    const twice = new Set<string>();
    for (const [i, k] of keys.entries()) {
      if (k === null)
        continue;
      if (rows.has(k))
        twice.add(k);
      rows.set(k, i);
    }
    return {keys, rows, twice};
  };
  const mine = rowsOf(ours);
  const theirs = rowsOf(server);
  const read = rowsOf(base);
  app.joining++;
  try {
    // a column we lack, or hold only as a layout's empty placeholder, comes as the server has it
    const added = new Set<string>();
    for (const col of server.columns) {
      const own = ours.col(col.name);
      if (own !== null && !isEmptyPlaceholder(own))
        continue;
      if (own === null)
        ours.columns.addNew(col.name, col.type);
      else if (own.type !== col.type)
        ours.changeColumnType(col.name, col.type as DG.COLUMN_TYPE);
      if (col.semType)
        ours.col(col.name)!.semType = col.semType;
      added.add(col.name);
    }
    for (let i = 0; i < server.rowCount; i++) {
      const k = theirs.keys[i];
      if (k === null)
        continue;
      const row = mine.rows.get(k);
      if (row === undefined) {
        ours.rows.addNew(null, true);
        for (const col of server.columns)
          setCorrectedValue(ours.col(col.name)!, col, ours.rowCount - 1, i);
        continue;
      }
      const baseRow = read.rows.get(k);
      const threeWay = baseRow !== undefined && !mine.twice.has(k) && !theirs.twice.has(k) && !read.twice.has(k);
      for (const col of server.columns) {
        if (col.name === molColName)
          continue;
        const target = ours.col(col.name)!;
        if (added.has(col.name)) {
          setCorrectedValue(target, col, row, i);
          continue;
        }
        const baseCol = base?.col(col.name);
        const before = baseCol ? baseCol.getString(baseRow!) : '';
        if (threeWay && target.getString(row) === before && target.getString(row) !== col.getString(i))
          setCorrectedValue(target, col, row, i);
      }
    }
  } finally {
    app.joining--;
  }
  ours.fireValuesChanged();
}

/**
 * Runs the implicit "auto-merge on campaign open" pass: if the campaign has a saved
 * mergeConfig (`cfg`) with `autoMergeOnOpen` set and a `filePath` from the file share, re-read the
 * file and merge it into the campaign the open reads (`target`). Compute is force-disabled for this pass
 * and the campaign is NOT saved — the merged state lives in memory until the user saves
 * it themselves. Failures are logged and surfaced as a warning balloon; they never bubble
 * up to the open flow so the campaign still opens.
 */
export async function autoMergeOnCampaignOpen<T extends HitDesignTemplate>(
  app: HitDesignApp<T>, cfg: HitDesignMergeConfig | undefined, target: MergeTarget,
): Promise<void> {
  if (!cfg?.autoMergeOnOpen || !cfg.filePath) return;
  try {
    const exists = await grok.dapi.files.exists(cfg.filePath);
    if (!exists) {
      grok.shell.warning(`Auto-merge: source file not found or not accessible: ${cfg.filePath}`);
      return;
    }
    const fi = DG.FileInfo.fromString(cfg.filePath, '');
    const df = await readFileInfoToDataFrame(fi);
    if (!df) {
      grok.shell.warning('Auto-merge: could not read the saved source file; merge skipped.');
      return;
    }
    await df.meta.detectSemanticTypes();
    await grok.data.detectSemanticTypes(df);
    // Compute pipeline never runs on the implicit pass — only structural data flows in.
    const passCfg: HitDesignMergeConfig = {...cfg, runComputeOnNewRows: false};
    await mergeIntoCampaign(app, df, passCfg, target);
  } catch (e) {
    _package.logger.error(e);
    grok.shell.warning('Auto-merge failed; opening the campaign without merging.');
  }
}
