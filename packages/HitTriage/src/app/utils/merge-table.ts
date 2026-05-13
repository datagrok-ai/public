/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {HitDesignApp} from '../hit-design-app';
import {HitDesignMergeConfig, HitDesignMergeMode, HitDesignTemplate} from '../types';
import {HitDesignMolColName, TileCategoriesColName, ViDColName, ViDSemType} from '../consts';
import {calculateCellValues} from './calculate-single-cell';
import {registerMolsBatch} from './molreg';
import {_package} from '../../package';

const MOL_NAME_HINTS = ['molecule', 'smiles', 'mol', 'structure', 'canonical_smiles', 'canonicalsmiles'];
const VID_NAME_HINTS = ['v-id', 'vid', 'v_id'];

export function guessMoleculeColumn(df: DG.DataFrame): DG.Column | null {
  const cols = df.columns.toList();
  const semCols = cols.filter((c) => c.semType === DG.SEMTYPE.MOLECULE);
  if (semCols.length > 0) {
    const named = semCols.find((c) => MOL_NAME_HINTS.includes(c.name.toLowerCase()));
    return named ?? semCols[0];
  }
  const named = cols.find((c) => MOL_NAME_HINTS.includes(c.name.toLowerCase()));
  if (named) return named;
  return cols.find((c) => c.type === DG.COLUMN_TYPE.STRING) ?? null;
}

export function guessVidColumn(df: DG.DataFrame): DG.Column | null {
  const cols = df.columns.toList();
  const named = cols.find((c) => VID_NAME_HINTS.includes(c.name.toLowerCase()));
  if (named) return named;
  const semCol = cols.find((c) => c.semType === ViDSemType);
  if (semCol) return semCol;
  return null;
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
  for (const col of incomingDf.columns.toList()) {
    // In SMILES mode, never carry over a V-iD-like column from the incoming table:
    // campaign V-iDs are authoritative.
    if (cfg.mode === 'smiles' && col.name === ViDColName) continue;

    let target = col.name;
    if (cfg.molColName && col.name === cfg.molColName) target = campMolColName;
    if (cfg.mode === 'vid' && cfg.vidColName && col.name === cfg.vidColName) target = ViDColName;

    sources.push({col, target, isComputed: false});
  }
  if (computedDf) {
    for (const col of computedDf.columns.toList()) {
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
 * Merges `incomingDf` into the current campaign of `app` based on the supplied config.
 * Supports both SMILES-based and V-iD-based joins, optional addition of new rows,
 * and optional execution of the campaign compute pipeline on the incoming molecules.
 */
export async function mergeIntoCampaign<T extends HitDesignTemplate>(
  app: HitDesignApp<T>,
  incomingDf: DG.DataFrame,
  cfg: HitDesignMergeConfig,
): Promise<void> {
  if (!app.dataFrame) throw new Error('Campaign has no data');
  const campDf = app.dataFrame;
  const campMolColName = app.molColName;

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
    const compute = app.campaign?.template?.compute ?? app.template?.compute;
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

  // Add columns missing in the campaign
  const newlyAdded = new Set<string>();
  for (const s of sources) {
    if (!campDf.columns.contains(s.target)) {
      const nc = campDf.columns.addNew(s.target, s.col.type);
      if (s.col.semType) nc.semType = s.col.semType;
      newlyAdded.add(s.target);
    }
  }

  app.isJoining = true;
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
    const newRowStartIdx = campDf.rowCount;
    if (unmatched.length > 0) {
      const hasStages = (app.stages?.length ?? 0) > 0 && campDf.columns.contains(TileCategoriesColName);
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
            stageCol.set(newCampRow, app.stages[0], false);
        }
      }

      // Re-register new molecules so they get fresh, dictionary-backed V-iDs.
      // This deliberately overwrites any V-iD value carried over from the incoming table
      // (incoming V-iDs are not authoritative for the campaign).
      if (cfg.molColName) {
        const incMolCol = incomingDf.col(cfg.molColName);
        const campVidCol = campDf.col(ViDColName);
        if (incMolCol && campVidCol && app.campaignId) {
          const newMols = unmatched.map((i) => incMolCol.get(i)).filter((m) => !!m && m !== '');
          if (newMols.length > 0) {
            try {
              const vidMap = await registerMolsBatch(newMols, app.campaignId, app.appName);
              for (let j = 0; j < unmatched.length; j++) {
                const mol = incMolCol.get(unmatched[j]);
                const vid = mol ? vidMap.get(mol) : undefined;
                if (vid) campVidCol.set(newRowStartIdx + j, vid, false);
              }
            } catch (e) {
              _package.logger.error(e);
            }
          }
        }
      }
    }
  } finally {
    setTimeout(() => {app.isJoining = false;}, 500);
  }

  campDf.fireValuesChanged();
}
