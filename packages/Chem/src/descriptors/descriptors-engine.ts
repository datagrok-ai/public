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
import {CruxWorkerApi, DescriptorValues} from '../crux/crux-api';

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
  const input = molecules.toList() as (string | null)[];
  const result = await (await getCruxService()).calcDescriptors(input, descriptors);
  if (unreadableForTests) {
    input.forEach((m, i) => {
      if (m && unreadableForTests!.has(m))
        Object.values(result).forEach((values) => values.value[i] = null);
    });
  }
  await fromContainerWhereUnread(input, descriptors, result);
  return Object.entries(result).map(([name, {type, value}]) => DG.Column.fromList(type, name, value));
}

/**
 * The rows with a molecule whose Crux values are all null (a molecule Crux's reader rejects) get the container's
 * values, where the container reads it; such a row is counted and logged as a fallback, so that a gap in Crux's reader
 * never shows as empty cells without a word. A molecule neither reads stays empty, and is no fallback.
 */
async function fromContainerWhereUnread(input: (string | null)[], descriptors: string[],
  result: {[name: string]: DescriptorValues}): Promise<void> {
  const values = Object.values(result);
  const rows = input.flatMap((m, i) => m && values.length && values.every((v) => v.value[i] === null) ? [i] : []);
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
