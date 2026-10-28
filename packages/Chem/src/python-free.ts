/**
 * Chem's server-side Python in Crux mode ("Chemical Engine": Crux; crux-core decision chem-python-free): the scripts
 * Curate, Mutate, Generate Conformers, Gasteiger Partial Charges (the panel), Butina clustering, InChI to molecule,
 * Calculate logP, Murcko Scaffolds and Filter by Catalogs, and the chemistry of Map Identifiers (the chem-chem
 * container's conversions) run in the browser on Chem's Crux file, each the script's or the container's own steps on
 * crux's ports, checked against them: the column work on the Crux workers (crux/jobs/python-free.ts), a panel's
 * one molecule on the main thread's module, InChIs by crux-sketch's InChI module.
 *
 * Every function branches once, by `byEngineAsync` (engine/chem-engine.ts): the script or the container is the RDKit
 * mode and the fallback, unchanged. The scripts are renamed with the suffix `Python` ({@link PYTHON_SCRIPTS}; the
 * Gasteiger script keeps its name and panel role, hidden behind the package's panel) and run without their outputs
 * processed: the package function that calls one returns its outputs, joined or shown as the script's were.
 */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {byEngineAsync} from './engine/chem-engine';
import {loadCrux} from './engine/crux-module';
import {getCruxService} from './crux/crux-searches';
import {cruxReadsCxExtension} from './crux/crux-service';
import {inchisFromInputs, mapYielding} from './crux/crux-inchi';
import {loadCruxSketch} from './crux/crux-sketcher';
import {BUTINA_FP_WORDS, ButinaNeighbours} from './crux/jobs/python-free';
import * as container from './docker/api';

/** The scripts' names: the package functions took theirs (`Curate`, `Mutate`, ...); the Gasteiger script keeps its
 * own (its graphics output has callers: a calculated column, Flow), only the panel is a package function. */
export const PYTHON_SCRIPTS = {curate: 'CuratePython', mutate: 'MutatePython', conformers: 'GenerateConformersPython',
  gasteiger: 'ChemistryGasteigerPartialCharges', butina: 'ButinaMoleculesClusteringPython',
  inchiToMol: 'InchiToMolPython', logP: 'CalculateLogPPython', murcko: 'MurckoScaffoldsPython',
  filters: 'FilterByCatalogsPython'};

/** filter_catalogs.py's catalogs (its `catalog` input's choices). */
export const FILTER_CATALOGS = ['BRENK', 'NIH', 'PAINS_A', 'PAINS_B', 'PAINS_C', 'ZINC'];

/** The Crux workers' time per molecule (a chunk's limit: 30 s plus this per molecule; logP: the descriptor jobs'), and
 * the one-molecule jobs' deadlines; past them the script runs. */
const CURATE_MS_PER_MOLECULE = 50;
const CURATE_TAUTOMERS_MS_PER_MOLECULE = 500;
const BUTINA_MS_PER_MOLECULE = 20;
/** Butina's neighbours of a row: a Tanimoto per row of the column, well under this each. */
const BUTINA_MS_PER_PAIR = 0.001;
const MURCKO_MS_PER_MOLECULE = 20;
const FILTERS_MS_PER_MOLECULE = 50;
const IDENTIFIERS_MS_PER_MOLECULE = 20;
const MUTATE_TIMEOUT_MS = 60000;
const CONFORMERS_TIMEOUT_MS = 60000;

/** This package's script `name` (by its exact name: never the package function of the same friendly name). */
export function pythonScript(name: string): DG.Func {
  const script = DG.Func.find({package: 'Chem', name}).find((f) => f.name === name);
  if (!script)
    throw new Error(`Chem: no script ${name}`);
  return script;
}

/** The script `name`'s first output, the script called with `params` without its outputs processed (no join, no
 * view). */
async function scriptOutput<T>(name: string, params: {[name: string]: any}): Promise<T> {
  return (await pythonScript(name).prepare(params).call(false, undefined, {processed: false}))
    .getOutputParamValue() as T;
}

/** Throws (the script or the container runs) where a molecule is a CXSMILES whose query features Crux's reader skips:
 * they read them. */
function assertCruxReads(molecules: (string | null)[]): void {
  const query = molecules.find((m) => !!m && !cruxReadsCxExtension(m));
  if (query !== undefined)
    throw new Error(`Crux does not read the query features of ${query}`);
}

/** Each distinct value of `col` once (its categories; '' for an empty cell), and each row's index in them, as the
 * column is now. */
function distinct(col: DG.Column): {values: string[], rows: ArrayLike<number>} {
  if (col.type === DG.TYPE.STRING)
    return {values: col.categories, rows: col.getRawData().slice(0, col.length)};
  const values = col.toList().map((v) => v === null || v === undefined ? '' : String(v));
  return {values, rows: Int32Array.from(values, (_, i) => i)};
}

/** A string column of molecules (SMILES unless `units`). */
function moleculeColumn(name: string, values: string[], units: string = DG.UNITS.Molecule.SMILES): DG.Column {
  const col = DG.Column.fromStrings(name, values);
  col.semType = DG.SEMTYPE.MOLECULE;
  col.meta.units = units;
  return col;
}

/** A job of crux/jobs/python-free*.ts over `items` (each chunk the job's first argument, then `args`) on the Crux
 * workers, the chunks' results joined. */
async function onWorkers<T>(op: string, items: (string | null)[], msPerItem: number, args: unknown[] = []):
  Promise<T[]> {
  const parts = await (await getCruxService()).jobInChunks<ArrayLike<T>>(op, items, msPerItem, args);
  return parts.flatMap((part) => Array.from(part));
}

// ---- Curate ----

/** curate.py's inputs. */
export interface CurateOptions {
  kekulization: boolean, normalization: boolean, reionization: boolean, neutralization: boolean,
  tautomerization: boolean, mainFragment: boolean,
}

/** curate.py's `curated` table: `curated_molecule`, a SMILES per row ('' where the molecule does not read). */
export function curate(data: DG.DataFrame, molecules: DG.Column, options: CurateOptions): Promise<DG.DataFrame> {
  return byEngineAsync('Curate', async () => {
    const {values, rows} = distinct(molecules);
    assertCruxReads(values);
    const ms = options.tautomerization ? CURATE_TAUTOMERS_MS_PER_MOLECULE : CURATE_MS_PER_MOLECULE;
    // the script's SMILES per molecule, null where it raises
    const curated = await onWorkers<string | null>('chemCurate', values, ms, [JSON.stringify(options)]);
    const raises = curated.findIndex((v) => v === null);
    if (raises >= 0) // the script's run fails on it: the script runs, for its error
      throw new Error(`chemCurate: the script raises on ${values[raises]}`);
    const col = moleculeColumn('curated_molecule', Array.from(rows, (k) => curated[k] as string));
    return DG.DataFrame.fromColumns([col]);
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.curate, {data, molecules, ...options}));
}

// ---- Mutate ----

/** A seed for the random mutations: the kernel's `random` is never seeded, so each call draws anew. */
const randomSeed = (): number => crypto.getRandomValues(new Uint32Array(1))[0];

/** mutate.py's `mutations` table, of one molecule column. */
export function mutate(molecule: string, steps: number, randomize: boolean, maxRandomResults: number):
  Promise<DG.DataFrame> {
  return byEngineAsync('Mutate', async () => {
    assertCruxReads([molecule]);
    const options = JSON.stringify({steps, randomize, maxRandomResults, seed: randomSeed()});
    const mutations = await (await getCruxService()).idleWorkerJob<string[]>('chemMutate', [molecule ?? '', options],
      MUTATE_TIMEOUT_MS);
    const df = DG.DataFrame.fromColumns([moleculeColumn('mutations', mutations)]);
    df.name = 'mutations';
    return df;
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.mutate, {molecule, steps, randomize, maxRandomResults}));
}

// ---- Generate Conformers ----

/** One row of generate_conformers.py's table, as chem-wasm's `chemConformers` gives it. */
interface ConformerRow {smiles: string, molblock: string, conformer: number, energy: number | null,
  rmsd: number | null}

/** generate_conformers.py's `conformers` table: `smiles`, `molblock` (with hydrogens), `conformer`, `energy`
 * (MMFF94), `rmsd` (to the first conformer). */
export function generateConformers(molecule: string, numConformers: number, optimize: boolean,
  rmsThreshold: number, maxAttempts: number, seed: number): Promise<DG.DataFrame> {
  return byEngineAsync('Generate Conformers', async () => {
    assertCruxReads([molecule]);
    const options = JSON.stringify({numConformers, optimize, rmsThreshold, maxAttempts, randomSeed: seed,
      timeoutMs: CONFORMERS_TIMEOUT_MS});
    const {rows} = await (await getCruxService()).idleWorkerJob<{rows: ConformerRow[]}>('chemConformers',
      [molecule ?? '', options], CONFORMERS_TIMEOUT_MS);
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('smiles', rows.map((r) => r.smiles)),
      moleculeColumn('molblock', rows.map((r) => r.molblock), DG.UNITS.Molecule.MOLBLOCK),
      DG.Column.fromInt32Array('conformer', Int32Array.from(rows, (r) => r.conformer)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'energy', rows.map((r) => r.energy)),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'rmsd', rows.map((r) => r.rmsd)),
    ]);
    df.name = 'conformers';
    return df;
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.conformers, {molecule, num_conformers: numConformers, optimize,
    rms_threshold: rmsThreshold, max_attempts: maxAttempts, random_seed: seed}));
}

// ---- Gasteiger Partial Charges ----

/** The script panel's picture, as the platform shows a script's graphics. */
function picture(url: string): HTMLElement {
  const div = ui.div([], 'grok-scripting-image-container-info-panel');
  div.style.backgroundImage = `url('${url}')`;
  return div;
}

/** A script's graphics output as the picture's URL, null where the script drew nothing. The platform hands the output
 * over as the element it shows the picture in (`div.grok-scripting-image-container`, the picture its background image),
 * not as the PNG's base64; a base64 string is taken too. */
function graphicsUrl(output: unknown): string | null {
  if (typeof output === 'string')
    return output ? `data:image/png;base64,${output}` : null;
  if (!(output instanceof HTMLElement))
    return null;
  const url = /url\(["']?(.*?)["']?\)/.exec(output.style.backgroundImage)?.[1];
  return url && !/,(null)?$/.test(url) ? url : null; // no picture: the platform's URL of a null value
}

/** gasteiger_charges.py's picture of `mol`'s Gasteiger charges (a similarity map, 400 x 400; a PNG, Crux's drawn from
 * its SVG), as an element; an empty one where the script draws nothing (the molecule does not read). */
export function gasteigerPicture(mol: string, contours: number): Promise<HTMLElement> {
  return byEngineAsync('Gasteiger Partial Charges', async () => {
    assertCruxReads([mol]);
    if (!Number.isInteger(contours) || contours < 0) // the script raises; chem-wasm's u32 would wrap
      throw new Error(`Gasteiger Partial Charges: contours ${contours}`);
    const crux = await loadCrux('Gasteiger Partial Charges');
    const map = JSON.parse(crux.chemGasteigerMap(mol, contours, 400, 400)) as {svg: string} | null;
    if (!map)
      return ui.div();
    return picture(await svgToPng(map.svg, 400, 400));
  }, async () => {
    const url = graphicsUrl(await scriptOutput<unknown>(PYTHON_SCRIPTS.gasteiger, {mol, contours}));
    return url ? picture(url) : ui.div();
  });
}

/** `svg` drawn on a canvas of its size, as a PNG data URL (the script's picture is a PNG; a contour map's SVG holds a
 * rectangle per grid cell, a megabyte or more of them). */
async function svgToPng(svg: string, width: number, height: number): Promise<string> {
  const url = URL.createObjectURL(new Blob([svg], {type: 'image/svg+xml'}));
  try {
    const img = new Image(width, height);
    await new Promise<void>((resolve, reject) => {
      img.onload = () => resolve();
      img.onerror = () => reject(new Error('Gasteiger Partial Charges: the picture did not load'));
      img.src = url;
    });
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    canvas.getContext('2d')!.drawImage(img, 0, 0, width, height);
    return canvas.toDataURL('image/png');
  } finally {
    URL.revokeObjectURL(url);
  }
}

/** The Gasteiger Partial Charges panel: the picture, the number of contour lines under the settings icon (the
 * script panel's parameter). */
export function gasteigerWidget(mol: string, contours: number): DG.Widget {
  const host = ui.div([ui.loader()]);
  let latest = 0;
  const draw = async (n: number) => {
    const call = ++latest;
    let shown: HTMLElement;
    try {
      shown = await gasteigerPicture(mol, n);
    } catch (e: any) {
      shown = ui.divText(String(e?.message ?? e), 'chem-gasteiger-error');
    }
    if (call !== latest) // a later count's picture is on its way
      return;
    ui.empty(host);
    host.append(shown);
  };
  const input = ui.input.int('Contours', {value: contours, onValueChanged: (n) => draw(n ?? contours)});
  const settings = ui.div([input.root]);
  settings.style.display = 'none';
  const icon = ui.iconFA('cog', () => settings.style.display = settings.style.display === 'none' ? '' : 'none',
    'Contour lines');
  draw(contours);
  return new DG.Widget(ui.divV([ui.divH([icon], {style: {justifyContent: 'flex-end'}}), settings, host]));
}

// ---- InChI to molecule ----

/** inchiToMol.py's `MolToSmiles(MolFromInchi(id))`: in Crux mode by the Crux sketcher's InChI module (crux-sketch's
 * `inchiToMol`: RDKit's InchiToMol on crux's port of the InChI library, the module Chem already ships). */
export function inchiToMol(id: string): Promise<string> {
  return byEngineAsync('inchiToMol', async () => (await loadCruxSketch()).inchiToMol(id),
    () => scriptOutput<string>(PYTHON_SCRIPTS.inchiToMol, {id}));
}

// ---- Butina clustering ----

/**
 * RDKit's `Butina.ClusterData` over the rows' neighbours (each row's, itself included, in index order): rows by
 * their neighbour count, most first, ties the later row first; each unassigned row of more than its own neighbour
 * starts a cluster of its unassigned neighbours; the rest are clusters of their own in that order. The cluster
 * number of each row, in the order made.
 */
export function butinaClusters(counts: Int32Array, neighbours: (row: number) => ArrayLike<number>): Int32Array {
  const n = counts.length;
  const order = Array.from({length: n}, (_, i) => i).sort((a, b) => counts[b] - counts[a] || b - a);
  const ids = new Int32Array(n);
  const seen = new Uint8Array(n);
  let next = 0;
  let k = 0;
  for (; k < n && counts[order[k]] > 1; k++) {
    const idx = order[k];
    if (seen[idx])
      continue;
    seen[idx] = 1;
    ids[idx] = next;
    const nbrs = neighbours(idx);
    for (let j = 0; j < nbrs.length; j++) {
      if (!seen[nbrs[j]]) {
        seen[nbrs[j]] = 1;
        ids[nbrs[j]] = next;
      }
    }
    next++;
  }
  for (; k < n; k++) {
    if (!seen[order[k]])
      ids[order[k]] = next++;
  }
  return ids;
}

/** butina_cluster.py's `clusters` table: `cluster (Butina)`, each row's cluster number (Morgan r2/1024 fingerprints,
 * Tanimoto distance at most `distanceCutoff`). */
export function butinaCluster(data: DG.DataFrame, molecules: DG.Column, distanceCutoff: number):
  Promise<DG.DataFrame> {
  return byEngineAsync('Butina Molecules Clustering', async () => {
    const service = await getCruxService();
    const values = molecules.toList().map((v) => v ?? '') as string[];
    assertCruxReads(values);
    const parts = await service.jobInChunks<Uint32Array>('butinaFingerprints', values, BUTINA_MS_PER_MOLECULE);
    const fps = new Uint32Array(values.length * BUTINA_FP_WORDS);
    let at = 0;
    for (const part of parts) {
      fps.set(part, at);
      at += part.length;
    }
    const rows = Array.from({length: values.length}, (_, i) => i);
    const found = await service.jobInChunks<ButinaNeighbours>('butinaNeighbours', rows,
      Math.max(BUTINA_MS_PER_MOLECULE, values.length * BUTINA_MS_PER_PAIR), [fps, distanceCutoff]);
    const counts = new Int32Array(values.length);
    const lists: Int32Array[] = new Array(values.length);
    let row = 0;
    for (const part of found) {
      let offset = 0;
      for (let k = 0; k < part.counts.length; k++, row++) {
        counts[row] = part.counts[k];
        lists[row] = part.neighbours.subarray(offset, offset + part.counts[k]);
        offset += part.counts[k];
      }
    }
    const ids = butinaClusters(counts, (i) => lists[i]);
    return DG.DataFrame.fromColumns([DG.Column.fromInt32Array('cluster (Butina)', ids)]);
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.butina, {data, molecules, distanceCutoff}));
}

// ---- Calculate logP ----

/** calc_logP.py's `result` table: `clogP`, RDKit's `Descriptors.MolLogP` (Wildman-Crippen) per row, empty where the
 * cell is blank or its molecule does not read. In Crux mode the descriptor service's MolLogP (the same as the
 * script's, bit for bit), the molecules read as the script reads them. */
export function calculateLogP(table: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
  return byEngineAsync('Calculate logP', async () => {
    const {values, rows} = distinct(molecules);
    assertCruxReads(values);
    // the script's: a blank cell is no molecule; a MOL block's `\n` (two characters) a line break
    const asked = values.map((m) => !m.trim() ? '' : m.includes('M  END') ? m.replaceAll('\\n', '\n') : m);
    const {MolLogP} = await (await getCruxService()).calcDescriptors(asked, ['MolLogP']);
    const logp = Array.from(rows, (k) => asked[k] && !Number.isNaN(MolLogP.value[k]) ? MolLogP.value[k] : null);
    return DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'clogP', logp)]);
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.logP, {table, molecules}));
}

// ---- Murcko Scaffolds ----

/** murcko_scaffold.py's `scaffolds` table: per row the Murcko scaffold's SMILES (`MurckoScaffoldSmiles`, without
 * stereo), '' for a molecule without rings, empty where the molecule does not read as a SMILES. */
export function murckoScaffolds(data: DG.DataFrame, smiles: DG.Column): Promise<DG.DataFrame> {
  return byEngineAsync('Murcko Scaffolds', async () => {
    const {values, rows} = distinct(smiles);
    assertCruxReads(values);
    const scaffolds = await onWorkers<string | null>('chemMurckoScaffolds', values, MURCKO_MS_PER_MOLECULE);
    return DG.DataFrame.fromColumns([moleculeColumn('scaffolds', Array.from(rows, (k) => scaffolds[k] ?? ''))]);
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.murcko, {data, smiles}));
}

// ---- Filter by Catalogs ----

/** filter_catalogs.py's `filter` table: one bool column named `catalog`, true where an entry of RDKit's FilterCatalog
 * `catalog` matches the row's molecule (read as a SMILES). */
export function filterByCatalogs(data: DG.DataFrame, smiles: DG.Column, catalog: string): Promise<DG.DataFrame> {
  return byEngineAsync('Filter by Catalogs', async () => {
    if (!FILTER_CATALOGS.includes(catalog)) // the script raises (a KeyError): it runs, for its error
      throw new Error(`Filter by Catalogs: no catalog ${catalog}`);
    const {values, rows} = distinct(smiles);
    assertCruxReads(values);
    const flags = await onWorkers<number>('chemFilterCatalog', values, FILTERS_MS_PER_MOLECULE, [catalog]);
    return DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.BOOL, catalog,
      Array.from(rows, (k) => flags[k] === 1))]);
  }, () => scriptOutput<DG.DataFrame>(PYTHON_SCRIPTS.filters, {data, smiles, catalog}));
}

// ---- Map Identifiers' chemistry (widgets/identifiers.ts): the chem-chem container's conversions ----

/** A cell's text as the container gets it (`ids.toList()`'s value): null and non-strings as null. */
const text = (v: unknown): string | null => typeof v === 'string' ? v : null;

/** The container's conversion `byContainer` of `values`, or in Crux mode `byCrux`; `molecules`: the values are
 * molecules (a CXSMILES whose query features Crux's reader skips: the container's), else InChIs. */
function convert<R>(values: unknown[], molecules: boolean, byCrux: (values: (string | null)[]) => Promise<R[]>,
  byContainer: () => Promise<R[]>): Promise<R[]> {
  return byEngineAsync('Map Identifiers', async () => {
    const texts = values.map(text);
    if (molecules)
      assertCruxReads(texts);
    return byCrux(texts);
  }, byContainer);
}

/** The container's `/chem/molecules_to_canonical`: RDKit's canonical SMILES per molecule, '' where it does not read. */
export function smilesToCanonical(smiles: string[]): Promise<string[]> {
  return convert(smiles, true,
    (molecules) => onWorkers<string>('chemCanonicalSmiles', molecules, IDENTIFIERS_MS_PER_MOLECULE),
    () => container.smilesToCanonical(smiles));
}

/** The InChIs (the keys when `keys`) of the molecules as the container writes them: InChI's input by the Crux
 * workers, the library by crux-sketch's InChI module (crux/crux-inchi.ts); '' where the molecule does not read. */
async function inchisOf(molecules: (string | null)[], keys: boolean): Promise<string[]> {
  const inputs = await onWorkers<Uint8Array | null>('chemInchiInputs', molecules, IDENTIFIERS_MS_PER_MOLECULE);
  return (await inchisFromInputs(inputs.map((input) => input ?? ''), keys)).map((s) => s ?? '');
}

/** The container's `/chem/molecules_to_inchi`: the standard InChI per molecule, '' where it does not read. */
export function smilesToInchi(smiles: string[]): Promise<string[]> {
  return convert(smiles, true, (molecules) => inchisOf(molecules, false),
    () => container.smilesToInchi(smiles));
}

/** The container's `/chem/molecules_to_inchi_key`: the standard InChIKey per molecule, '' where it does not read. */
export function smilesToInchiKey(smiles: string[]): Promise<string[]> {
  return convert(smiles, true, (molecules) => inchisOf(molecules, true),
    () => container.smilesToInchiKey(smiles));
}

/** Whether `e` is crux-sketch's refusal of an InChI (no structure, an InChIKey, a structure RDKit cannot sanitize), not
 * its module's failure to load or run (`code`, or not a `SketchError`): the caller's container runs then. */
const refusedInchi = (e: any): boolean => !!e && typeof e === 'object' && e.format === 'inchi' && !e.code &&
  !(e instanceof Error);

/** `one` of each InChI by crux-sketch's InChI module (`none` for an empty one), yielding to the page now and then. */
function eachInchi<R>(inchis: (string | null)[], one: (inchi: string) => Promise<R>, none: R): Promise<R[]> {
  return mapYielding(inchis, (inchi) => inchi ? one(inchi) : none);
}

/** The container's `/chem/inchi_to_smiles`: `MolToSmiles(MolFromInchi(inchi))` per InChI, '' where RDKit makes no
 * molecule of it. In Crux mode crux-sketch's `inchiToMol` (RDKit's MolFromInchi on crux's InChI library; the SMILES
 * written under the new stereo perception, as Chem's: the container's legacy one writes a few otherwise). */
export function inchiToSmiles(inchi: string[]): Promise<string[]> {
  return convert(inchi, false, async (inchis) => {
    const sketch = await loadCruxSketch();
    return eachInchi(inchis, (s) => sketch.inchiToMol(s).catch((e) => {
      if (refusedInchi(e))
        return '';
      throw e;
    }), '');
  }, () => container.inchiToSmiles(inchi));
}

/** The container's `/chem/inchi_to_inchi_key`: the InChIKey per InChI, null for '' and where the library refuses it. */
export function inchiToInchiKey(inchi: string[]): Promise<string[]> {
  return convert(inchi, false, async (inchis) => {
    const sketch = await loadCruxSketch();
    return eachInchi(inchis, async (s) => (await sketch.inchiKeyFromInchi(s)) || null, null) as Promise<string[]>;
  }, () => container.inchiToInchiKey(inchi));
}
