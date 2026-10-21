/**
 * Descriptors by the "Chemical engine" setting (src/engine/README.md): in Crux mode the dialog's tree comes from the
 * Crux module and the values from the Crux workers, else (and whenever Crux fails) both from the chem-chem container.
 * Every descriptor entry point (the menu, its editor and transform, the vector function, the panel, `chemDescriptors`,
 * `chemDescriptorsTree`) goes through these two functions, with the container's results' shape.
 */
import * as DG from 'datagrok-api/dg';
import * as container from '../docker/api';
import {isCruxRequested} from '../engine/chem-engine';
import {cruxOrRdKitAsync, warnCruxFallback} from '../engine/fallback';
import {getCruxService} from '../crux/crux-searches';
import {loadCrux} from '../engine/crux-module';
import {CruxWorkerApi, DescriptorTiming, INT_NONE, PackedDescriptorValues, isNoneValue} from '../crux/crux-api';

let cruxTree: object | null = null;

/** The descriptors by group, as the container's `GET /chem/descriptors/tree` gives them. */
export async function getDescriptorsTree(): Promise<object> {
  if (!isCruxRequested())
    return container.getDescriptorsTree();
  // the tree is the module's static data: the main thread's instance gives it without starting the workers
  return cruxOrRdKitAsync('descriptorTree',
    async () => cruxTree ??= JSON.parse((await loadCrux('descriptors') as unknown as CruxWorkerApi).descriptorTree()),
    () => container.getDescriptorsTree());
}

/** Tests: the next tree in Crux mode comes from the loader again. */
export function forgetCruxTreeForTests(): void {
  cruxTree = null;
}

/** A column per descriptor (group names expand to their members), named, typed and with nulls as the container's. */
export async function calculateDescriptors(molecules: DG.Column, descriptors: string[]): Promise<DG.Column[]> {
  if (!isCruxRequested())
    return container.calculateDescriptors(molecules, descriptors);
  return cruxOrRdKitAsync('calcDescriptors',
    () => cruxDescriptors(molecules, descriptors),
    () => container.calculateDescriptors(molecules, descriptors));
}

async function cruxDescriptors(molecules: DG.Column, descriptors: string[]): Promise<DG.Column[]> {
  // the column's values as the container gets them: Crux reads them as chem.py does (a MOL block if it holds
  // "M  END", else SMILES; '' an empty molecule; null or unreadable: nulls)
  const timing: DescriptorTiming | null = timingOn ?
    {rows: molecules.length, names: descriptors.length, workers: 0, phases: {}, chunks: []} : null;
  let phase = performance.now();
  const lap = (name: string) => {
    if (timing) {
      const now = performance.now();
      timing.phases[name] = now - phase;
      phase = now;
    }
  };
  // each distinct molecule once (a string column's categories), its values then spread over its rows
  const byCategory = molecules.type === DG.COLUMN_TYPE.STRING;
  // (an empty cell is null, as toList() gives it, not the empty category '': Crux reads '' as an empty molecule)
  const input = (byCategory ? molecules.categories.map((c) => c === '' ? null : c) : molecules.toList()) as
    (string | null)[];
  const rowCategories = byCategory ? molecules.getRawData() as Int32Array : null;
  lap('read column');
  const service = await getCruxService();
  lap('pool');
  const result = await service.calcDescriptors(input, descriptors, timing ?? undefined);
  phase = performance.now();
  if (unreadableForTests) {
    input.forEach((m, i) => {
      if (m && unreadableForTests!.has(m))
        Object.values(result).forEach((values) => values.value[i] = values.type === 'int' ? INT_NONE : NaN);
    });
  }
  await fromContainerWhereUnread(input, descriptors, result);
  lap('unread rows');
  const columns = Object.entries(result).map(([name, values]) => descriptorColumn(name, values, rowCategories));
  lap('columns');
  if (timing) {
    lastTiming = timing;
    console.info(`Chem | Crux descriptors: ${timing.rows} rows, ${timing.names} names, ${timing.workers} workers; ` +
      Object.entries(timing.phases).map(([k, v]) => `${k} ${v.toFixed(0)} ms`).join(', '));
  }
  return columns;
}

/**
 * A descriptor's column, as DG.Column.fromList(type, name, values with nulls) makes it (its storage too), from the
 * typed values, spread over the rows by `rowCategories` (a row's category) when given.
 */
function descriptorColumn(name: string, values: PackedDescriptorValues, rowCategories: Int32Array | null): DG.Column {
  const length = rowCategories ? rowCategories.length : values.value.length;
  const at = rowCategories ? (row: number) => rowCategories[row] : (row: number) => row;
  if (values.type === 'int') {
    const array = new Int32Array(length);
    for (let row = 0; row < length; row++)
      array[row] = values.value[at(row)]; // INT_NONE is Datagrok's int null
    return DG.Column.fromInt32Array(name, array);
  }
  doubleStorage ??= DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, '', [0]).getRawData() instanceof Float64Array ?
    'float64' : 'float32';
  const array = doubleStorage === 'float64' ? new Float64Array(length) : new Float32Array(length);
  for (let row = 0; row < length; row++) {
    const v = values.value[at(row)];
    array[row] = Number.isNaN(v) ? DG.FLOAT_NULL : v;
  }
  return array instanceof Float64Array ? DG.Column.fromFloat64Array(name, array) :
    DG.Column.fromFloat32Array(name, array);
}
/** How DG.Column.fromList stores doubles (the columns made here match it). */
let doubleStorage: 'float32' | 'float64' | null = null;

let timingOn = false;
let lastTiming: DescriptorTiming | null = null;

/** Diagnosis: Crux descriptors calls are timed per phase and per chunk (and logged) while on. */
export function setDescriptorTimingForTests(on: boolean): void {
  timingOn = on;
  lastTiming = null;
}

/** The last timed call's phases and chunks (null when none was timed). */
export function lastDescriptorTiming(): DescriptorTiming | null {
  return lastTiming;
}

/**
 * The rows with a molecule whose Crux values are all null (a molecule Crux's reader rejects) get the container's
 * values, where the container reads it; such a row is counted and logged as a fallback, so that a gap in Crux's reader
 * never shows as empty cells without a word. A molecule neither reads stays empty, and is no fallback.
 */
async function fromContainerWhereUnread(input: (string | null)[], descriptors: string[],
  result: {[name: string]: PackedDescriptorValues}): Promise<void> {
  const values = Object.values(result);
  const rows = input.flatMap((m, i) => m && values.length && values.every((v) => isNoneValue(v, i)) ? [i] : []);
  if (!rows.length)
    return;
  let columns: DG.Column[];
  try {
    columns = await container.calculateDescriptors(DG.Column.fromStrings('molecules', rows.map((i) => input[i]!)),
      descriptors);
  } catch (e: any) {
    if (!containerWarned) {
      containerWarned = true;
      console.warn(`Chem | Crux could not read ${rows.length} molecules and the container is not reachable: their ` +
        `descriptors stay empty (e.g. ${input[rows[0]]?.slice(0, 80)}): ${e?.message ?? e}`);
    }
    return;
  }
  const read = new Set<number>();
  for (const column of columns) {
    const target = result[column.name];
    if (!target)
      continue;
    rows.forEach((row, k) => {
      const v = column.isNone(k) ? null : column.get(k);
      // (a container image before GROK-20505 gives some values as text)
      const n = v === null || v === '' ? null : typeof v === 'string' ? Number(v) : v;
      if (n !== null && !Number.isNaN(n)) {
        target.value[row] = n;
        read.add(row);
      }
    });
  }
  if (read.size) {
    warnCruxFallback('calcDescriptors', `Crux could not read ${read.size} molecules, the container calculated ` +
      `them, e.g. ${input[[...read][0]]?.slice(0, 80)}`);
  }
}
let containerWarned = false;

/** Tests: Crux's values for these molecules come back all null, as for a molecule its reader rejects (null: none). */
export function setCruxDescriptorGapForTests(molecules: string[] | null): void {
  unreadableForTests = molecules ? new Set(molecules) : null;
}
let unreadableForTests: Set<string> | null = null;
