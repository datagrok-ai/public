/**
 * Column jobs by the "Chemical engine" setting (README.md beside this file) that ran on RDKit in both modes, on the
 * RDKit pool or on the main thread: in Crux mode the Crux workers run them in chunks (crux/jobs/molecule-jobs.ts), the
 * RDKit path does the rows Crux leaves (a molecule it does not read, a call it refuses, a CXSMILES whose query features
 * its reader skips) and the whole call when Crux fails; in RDKit mode, the RDKit path as before.
 */
import * as DG from 'datagrok-api/dg';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {getCruxService} from '../crux/crux-searches';
import type {ElementCounts, MoleculeJobResult} from '../crux/jobs/molecule-jobs';
import {byEngineAsync} from './chem-engine';
import {rdkitForRows} from './fallback';
import {cruxInput} from './worker-jobs';

/** Per molecule, the Crux workers' time allowance of a job (crux-service.ts `inChunks`): CoordGen lays out a
 * macrocycle in seconds. */
const MS_PER_MOLECULE = {coordGen: 200, flatten: 20, beautify: 100, elements: 20, cores: 50, atomCounts: 5};

/** The molecule job `op` (crux/jobs/molecule-jobs.ts) over `molecules` on the Crux workers: each chunk's result and
 * the rows the RDKit side does (the chunks' `rest` and the skipped CXSMILES), as column rows. */
async function moleculeJob<R extends {rest: number[]}>(op: string, molecules: string[], msPerMolecule: number):
  Promise<{parts: R[], starts: number[], rest: number[]}> {
  const {input, skipped} = cruxInput(molecules);
  const starts: number[] = [];
  const parts = await (await getCruxService()).jobInChunks<R>(op, input, msPerMolecule, [], starts);
  return {parts, starts, rest: [...parts.flatMap((part, k) => part.rest.map((row) => starts[k] + row)), ...skipped]};
}

/**
 * `op`'s results over `molecules` on the Crux workers in Crux mode, `rdkit` (the RDKit path, over some molecules) for
 * the rows the workers leave and for the whole call when they fail; `rdkit` alone in RDKit mode. A row `rdkit` gives
 * `done` for is logged as a fallback.
 */
async function resultsByEngine<T>(call: string, op: string, molecules: string[], msPerMolecule: number,
  rdkit: (molecules: string[]) => Promise<T[]>, done: (result: T) => boolean): Promise<T[]> {
  return byEngineAsync(call, async () => {
    const {parts, rest} = await moleculeJob<MoleculeJobResult<T>>(op, molecules, msPerMolecule);
    const results = parts.flatMap((part) => part.results);
    await rdkitForRows(call, rest, (rows) => rdkit(rows.map((i) => molecules[i])),
      (patch, k, row) => done(results[row] = patch[k]), (row) => molecules[row]);
    return results;
  }, () => rdkit(molecules));
}

/** Transform | Recalculate Coordinates by CoordGen: `RdKitService.getCoordGenCoords(molecules)` by the engine. */
export function coordGenCoordsByEngine(molecules: string[]): Promise<string[]> {
  return resultsByEngine('CoordGen coordinates', 'coordGenCoords', molecules, MS_PER_MOLECULE.coordGen,
    async (rows) => (await getRdKitService()).getCoordGenCoords(rows), (r) => !!r);
}

/** Transform | Flatten Molecules: `RdKitService.flattenMolecules(molecules)` (SMILES without stereo) by the engine. */
export function flattenMoleculesByEngine(molecules: string[]): Promise<string[]> {
  return resultsByEngine('flattenMolecules', 'flattenMolecules', molecules, MS_PER_MOLECULE.flatten,
    async (rows) => (await getRdKitService()).flattenMolecules(rows), (r) => !!r);
}

/** `beautifyMols` (Bio's To Atomic Level): `RdKitService.beautifyMolsV3K(molecules)` by the engine. */
export function beautifyMolsByEngine(molecules: string[]): Promise<string[]> {
  return resultsByEngine('beautifyMols', 'beautifyMols', molecules, MS_PER_MOLECULE.beautify,
    async (rows) => (await getRdKitService()).beautifyMolsV3K(rows), (r) => !!r);
}

/**
 * The R-group dialog's Core column in Crux mode (analysis/r-group-analysis.ts): each core's MOL block as the dialog
 * writes it (null: none), from the Crux workers; `undefined` for a core the dialog's own RDKit code does (Crux reads
 * none from it). Null in RDKit mode, or when the workers fail (the dialog's code does every core).
 */
export async function coreMolblocksByCrux(cores: string[]): Promise<(string | null | undefined)[] | null> {
  return byEngineAsync('R-group cores', async () => {
    const {parts, rest} = await moleculeJob<MoleculeJobResult<string | null>>('coreMolblocks', cores,
      MS_PER_MOLECULE.cores);
    const results: (string | null | undefined)[] = parts.flatMap((part) => part.results.map((r) => r || null));
    for (const row of rest)
      results[row] = undefined;
    return results;
  }, async () => null);
}

/** Heavy-atom counts of structures (the SAR matrix's `cachedAtomCount`) by the engine: `rdkit` per structure in RDKit
 * mode (its one caller, the SAR matrix's prefetch, asks in Crux mode only), and for those Crux reads none from. */
export function atomCountsByEngine(structures: string[], rdkit: (s: string) => number): Promise<number[]> {
  return resultsByEngine('atom counts', 'atomCounts', structures, MS_PER_MOLECULE.atomCounts,
    async (rows) => rows.map(rdkit), (n) => n > 0);
}

/**
 * Elemental analysis (`getAtomsColumn` of utils/elemental-analysis-utils.ts, passed in as `rdkit`) by the engine: in
 * Crux mode the Crux workers count each molecule's atoms, `rdkit` the rows they leave (on the main thread, read as
 * the column's first row tells: SMILES, V3000 or V2000), the whole column when they fail; in RDKit mode, `rdkit`.
 */
export async function atomsColumnByEngine(molCol: DG.Column,
  rdkit: (molCol: DG.Column) => [Map<string, Int32Array>, number[]]): Promise<[Map<string, Int32Array>, number[]]> {
  return byEngineAsync('elemental analysis', async () => {
    const molecules = molCol.toList();
    const n = molecules.length;
    const {parts, starts, rest} = await moleculeJob<ElementCounts & {rest: number[]}>('elementalAnalysis', molecules,
      MS_PER_MOLECULE.elements);
    const charge = new Int32Array(n);
    const elements = new Map<string, Int32Array>([['Molecule Charge', charge]]);
    const elementColumn = (symbol: string): Int32Array => {
      if (!elements.has(symbol))
        elements.set(symbol, new Int32Array(n));
      return elements.get(symbol)!;
    };
    parts.forEach((part, k) => {
      for (const [symbol, counts] of Object.entries(part.elements))
        elementColumn(symbol).set(counts, starts[k]);
      charge.set(part.charges, starts[k]);
    });
    const invalid = new Array<number>(n);
    // the rows on RDKit's path, after the column's first row: it picks how the path reads them
    await rdkitForRows('elemental analysis', rest, async (rows) => {
      const restCol = DG.Column.fromStrings('molecules', [molecules[0], ...rows.map((i) => molecules[i])]);
      if (molCol.meta.units)
        restCol.meta.units = molCol.meta.units;
      return rdkit(restCol);
    }, ([restElements, restInvalid], k, row) => {
      for (const [symbol, counts] of restElements)
        elementColumn(symbol)[row] = counts[k + 1];
      if (restInvalid[k + 1] !== undefined) {
        invalid[row] = row;
        return false;
      }
      return !!molecules[row];
    }, (row) => molecules[row]);
    return [elements, invalid];
  }, async () => rdkit(molCol));
}
