import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {SEMTYPE} from '../utils/proteomics-types';
import {RunMeta, setRunMeta, getRunMeta} from './spc';

/** Describes which intensity columns belong to each experimental group. */
export interface GroupAssignment {
  group1: {name: string; columns: string[]};
  group2: {name: string; columns: string[]};
}

/** Store group assignments as a DataFrame tag. */
export function setGroups(df: DG.DataFrame, groups: GroupAssignment): void {
  df.setTag('proteomics.groups', JSON.stringify(groups));
}

/** Retrieve group assignments from DataFrame tag. Returns null if not set. */
export function getGroups(df: DG.DataFrame): GroupAssignment | null {
  const raw = df.getTag('proteomics.groups');
  if (!raw) return null;
  try {
    return JSON.parse(raw) as GroupAssignment;
  } catch {
    return null;
  }
}

/** Seed values for the Annotate Experiment dialog. */
export interface AnnotationSeed {
  group1Name: string;
  group1Cols: DG.Column[];
  group2Name: string;
  group2Cols: DG.Column[];
  /** Phase 16 D-01: persisted operator value, else parser-side seed, else undefined. */
  instrumentId?: string;
  /** Phase 16 D-01: ISO-8601 (date-only acceptable per RESEARCH gotcha #1). */
  acquisitionDatetime?: string;
}

/**
 * Seeds Annotate Experiment dialog inputs from existing `proteomics.groups` if present,
 * else returns the legacy Control/Treatment defaults. Stale column refs (no longer in the
 * DataFrame or in the `available` list) are silently dropped — protects parser-populated
 * assignments (e.g. Spectronaut DMD/WT) from being blanked when the user OKs the dialog
 * without re-selecting columns.
 */
export function seedAnnotationDialogInputs(df: DG.DataFrame, available: string[]): AnnotationSeed {
  const runMetaSeed = readRunMetaSeed(df);
  const existing = getGroups(df);
  if (!existing) {
    return {
      group1Name: 'Control', group1Cols: [],
      group2Name: 'Treatment', group2Cols: [],
      instrumentId: runMetaSeed?.instrument_id,
      acquisitionDatetime: runMetaSeed?.acquisition_datetime,
    };
  }
  const availableSet = new Set(available);
  const resolve = (names: string[]): DG.Column[] => names
    .filter((n) => availableSet.has(n))
    .map((n) => df.col(n))
    .filter((c): c is DG.Column => c !== null);
  return {
    group1Name: existing.group1.name || 'Control',
    group1Cols: resolve(existing.group1.columns),
    group2Name: existing.group2.name || 'Treatment',
    group2Cols: resolve(existing.group2.columns),
    instrumentId: runMetaSeed?.instrument_id,
    acquisitionDatetime: runMetaSeed?.acquisition_datetime,
  };
}

/** Reads the Phase 16 run-identity seed in priority order:
 *  1. Persisted `proteomics.spc_run_meta` (operator value from a prior dialog OK)
 *  2. Parser-side `proteomics.spc_run_meta_seed` (set by the Spectronaut PG parser)
 *  3. undefined (operator types both fields manually)
 */
function readRunMetaSeed(df: DG.DataFrame): RunMeta | undefined {
  const persisted = getRunMeta(df);
  if (persisted) return persisted;
  const raw = df.getTag('proteomics.spc_run_meta_seed');
  if (!raw) return undefined;
  try {
    return JSON.parse(raw) as RunMeta;
  } catch {
    return undefined;
  }
}

/** Shows a dialog for assigning intensity columns to two experimental groups. */
export function showAnnotationDialog(df: DG.DataFrame): void {
  const intensityColNames = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map((c) => c.name);

  const seed = seedAnnotationDialogInputs(df, intensityColNames);

  const group1Name = ui.input.string('Group 1 Name', {value: seed.group1Name});
  const group1Cols = ui.input.columns('Group 1', {table: df, available: intensityColNames, value: seed.group1Cols});
  const group2Name = ui.input.string('Group 2 Name', {value: seed.group2Name});
  const group2Cols = ui.input.columns('Group 2', {table: df, available: intensityColNames, value: seed.group2Cols});
  // Phase 16 D-01: two new optional run-identity inputs. Tooltips come VERBATIM from UI-SPEC.
  const instrumentInput = ui.input.string('Instrument', {
    value: seed.instrumentId ?? '',
    tooltipText: 'Free-text identifier for the instrument that ran this analysis (e.g. QExactive-01). ' +
      'Used to group runs in the SPC dashboard.',
  });
  const datetimeInput = ui.input.date('Acquisition datetime', {
    value: seed.acquisitionDatetime ? dayjs(seed.acquisitionDatetime) : null,
    tooltipText: 'When the samples were acquired on the instrument. ' +
      'Used to order runs on the SPC trend chart — not the file-import time.',
  } as any);

  ui.dialog('Annotate Experiment')
    .add(group1Name)
    .add(group1Cols)
    .add(group2Name)
    .add(group2Cols)
    .add(instrumentInput)
    .add(datetimeInput)
    .onOK(() => {
      const g1: DG.Column[] = group1Cols.value ?? [];
      const g2: DG.Column[] = group2Cols.value ?? [];

      if (g1.length === 0 || g2.length === 0) {
        grok.shell.warning('Select at least one column per group');
        return;
      }

      const groups: GroupAssignment = {
        group1: {name: group1Name.value, columns: g1.map((c) => c.name)},
        group2: {name: group2Name.value, columns: g2.map((c) => c.name)},
      };
      const instrumentValue = (instrumentInput.value ?? '').trim();
      const datetimeValue: any = datetimeInput.value;
      const isoDatetime = datetimeValue && typeof datetimeValue.toISOString === 'function'
        ? datetimeValue.toISOString()
        : (datetimeValue instanceof Date ? datetimeValue.toISOString() : '');
      applyAnnotation(df, {
        group1: groups.group1,
        group2: groups.group2,
        runMeta: (instrumentValue.length > 0 || isoDatetime.length > 0)
          ? {instrument_id: instrumentValue, acquisition_datetime: isoDatetime}
          : undefined,
      });
      grok.shell.info(`Groups assigned: ${g1.length} + ${g2.length} samples`);
    })
    .show();
}

/** Programmatic Annotate Experiment OK-path. Used by Plan 16-01's
 *  `SPC:annotation_dialog_persists_run_meta` test and by the dialog above to
 *  share the persist-on-OK flow. setGroups + setRunMeta in lockstep so the
 *  test never has to instantiate a real dialog. */
export function applyAnnotation(
  df: DG.DataFrame,
  payload: {
    group1: GroupAssignment['group1'];
    group2: GroupAssignment['group2'];
    runMeta?: RunMeta;
  },
): void {
  setGroups(df, {group1: payload.group1, group2: payload.group2});
  if (payload.runMeta) setRunMeta(df, payload.runMeta);
}
