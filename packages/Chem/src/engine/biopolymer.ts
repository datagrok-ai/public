/**
 * Chem's macromolecule assembler (`Chem:getMacromoleculeAssembler`, chem-meta's IMacromoleculeAssembler; Bio-on-Crux
 * PLAN §3.3, §5.3): Bio's parsed sequences assembled from the templates of a monomer library by Crux, in both chemical
 * engines, as `molzip` is (Crux loads in RDKit mode too; there is no RDKit path, so a call fails when Crux does not
 * load). A column (`assemble`, a long list for `markStereoAbsolute`) runs on the Crux workers (CruxService); one row
 * (`assembleOne`, `toV2000`, a few MOL blocks) on the main thread's module, which compiles each library it is asked for
 * once, as each worker does. The synchronous twins (`assembleOneSync`, `assembleSync`, `markStereoAbsoluteSync`,
 * `toV2000Sync`) serve Bio's legacy synchronous APIs on that module: they need `ready()` first (the module loaded,
 * the library registered), and compile a library on first use there.
 */
import {AsmBatch, AsmLayout, AsmSequence, AssembleOneResult, AssembleOptions, AssembleResult, IMacromoleculeAssembler}
  from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {getCruxService} from '../crux/crux-searches';
import {CruxAssembledBatch, CruxAssembledRow, CruxBiopolymerApi, CruxMonomerLibrary, assembleOptionsJson,
  cruxBiopolymer} from '../crux/crux-api';
import {cruxLoadFailed, loadCrux, noteCruxError} from './crux-module';

/** The libraries kept by key, the least recently used dropped first (Bio keeps a few, each a new key per change). */
const MAX_LIBRARIES = 8;
/** The libraries compiled on the main thread, the least recently used freed first. */
const MAX_COMPILED = 4;
/** Up to this many MOL blocks `markStereoAbsolute` does on the main thread; more on the Crux workers. */
const MAIN_THREAD_MOLBLOCKS = 100;
/** Counts the worker keys given out: unique in the session, whichever assembler gives them (the pool is shared). */
let lastWorkerKey = 0;

/** `call()` on the main thread's module; a trap there is noted (Crux is off for the session: its instance may be
 * broken) and thrown. */
function onMainThread<T>(call: () => T): T {
  try {
    return call();
  } catch (e) {
    noteCruxError(e);
    throw e instanceof Error ? e : new Error(String(e)); // the binding throws strings for refusals
  }
}

/** A copy of `positions` as an Int32Array of its own (not a view of the wasm memory). */
const int32s = (positions: Int32Array | number[]): Int32Array => Int32Array.from(positions);

export class CruxMacromoleculeAssembler implements IMacromoleculeAssembler {
  // By the caller's key, the most recently used last; `id` is the workers' key, a new one for each JSON, so a key set
  // again with another library is compiled again.
  private libraries = new Map<string, {json: string, id: string}>();
  // The main thread's compiled libraries by worker key, the most recently used last, and the module they are of.
  private compiled = new Map<string, CruxMonomerLibrary>();
  private compiledOf: unknown = null;
  // The main thread's module once loaded, for the synchronous twins (set by every async main-thread call and ready())
  private api: CruxBiopolymerApi | null = null;

  async setLibrary(key: string, json: string): Promise<void> {
    if (typeof json !== 'string' || json === '')
      throw new Error(`Chem | the monomer library ${key} is not Pistoia JSON text`);
    const known = this.libraries.get(key);
    this.libraries.delete(key);
    if (known?.json === json) {
      this.libraries.set(key, known);
      return;
    }
    if (known)
      this.freeCompiled(known.id);
    // a worker key of Chem's own: a caller's key could hold a word the pool reads as a trap's (CruxService isTrap)
    this.libraries.set(key, {json, id: `library ${++lastWorkerKey}`});
    while (this.libraries.size > MAX_LIBRARIES) {
      const [oldest, entry] = this.libraries.entries().next().value!;
      this.freeCompiled(entry.id);
      this.libraries.delete(oldest);
    }
  }

  async assemble(batch: AsmBatch, libraryKey: string, options: AssembleOptions): Promise<AssembleResult> {
    const {json, id} = this.library(libraryKey);
    return (await getCruxService()).assembleMacromolecules(batch, id, json, options);
  }

  async assembleOne(seq: AsmSequence, libraryKey: string, options: AssembleOptions): Promise<AssembleOneResult> {
    const {json, id} = this.library(libraryKey);
    const library = await this.mainThreadLibrary(id, json);
    const options1 = assembleOptionsJson(options, true);
    return assembledOne(onMainThread(() => library.assembleOne(JSON.stringify(seq), options1)), options);
  }

  /** An empty cell (null, '', a hole of `toList()`) stays as it is ('' stays '', the rest null). */
  async markStereoAbsolute(molblocks: (string | null)[]): Promise<(string | null)[]> {
    const {result, rows, sent} = stereoInput(molblocks);
    if (sent.length === 0)
      return result;
    let marked: (string | null)[];
    if (sent.length <= MAIN_THREAD_MOLBLOCKS) {
      const crux = await this.mainThreadApi();
      marked = onMainThread(() => crux.chemMarkStereoAbsolute(sent));
    } else
      marked = await (await getCruxService()).markStereoAbsolute(sent);
    rows.forEach((row, k) => result[row] = marked[k] ?? null);
    return result;
  }

  async toV2000(molblock: string): Promise<string> {
    const crux = await this.mainThreadApi();
    return onMainThread(() => crux.chemToV2000(molblock));
  }

  // -- the synchronous twins (main thread only) --

  async ready(libraryKey: string, json: string): Promise<void> {
    await this.setLibrary(libraryKey, json);
    const {json: libJson, id} = this.library(libraryKey);
    await this.mainThreadLibrary(id, libJson);
  }

  assembleOneSync(seq: AsmSequence, libraryKey: string, options: AssembleOptions): AssembleOneResult {
    const {json, id} = this.library(libraryKey);
    const library = this.mainThreadLibrarySync(id, json);
    const options1 = assembleOptionsJson(options, true);
    return assembledOne(onMainThread(() => library.assembleOne(JSON.stringify(seq), options1)), options);
  }

  assembleSync(batch: AsmBatch, libraryKey: string, options: AssembleOptions): AssembleResult {
    const {json, id} = this.library(libraryKey);
    const library = this.mainThreadLibrarySync(id, json);
    const n = batch.rowStart.length - 1;
    const r: CruxAssembledBatch = onMainThread(() => library.assemble(JSON.stringify(batch.symbols), batch.rowStart,
      batch.monSymbol, batch.monRole, batch.monPosition, batch.chainStart, batch.rowChainStart, batch.conn,
      batch.rowConnStart, Array.from(batch.rowError, (e) => e ?? null), assembleOptionsJson(options, false)));
    const errors = Array.from({length: n}, (_, i) => r.errors[i] ?? null);
    const molblocks = Array.from({length: n}, (_, i) => errors[i] ? null : r.molblocks[i] ?? null);
    return {molblocks, layouts: Uint8Array.from(r.layouts), errors}; // a copy: not a view of the wasm memory
  }

  markStereoAbsoluteSync(molblocks: (string | null)[]): (string | null)[] {
    const {result, rows, sent} = stereoInput(molblocks);
    if (sent.length === 0)
      return result;
    const crux = this.syncApi();
    const marked = onMainThread(() => crux.chemMarkStereoAbsolute(sent));
    rows.forEach((row, k) => result[row] = marked[k] ?? null);
    return result;
  }

  toV2000Sync(molblock: string): string {
    const crux = this.syncApi();
    return onMainThread(() => crux.chemToV2000(molblock));
  }

  /** Tests: the key the workers know library `key` by. */
  workerKeyForTests(key: string): string | undefined {
    return this.libraries.get(key)?.id;
  }

  private library(key: string): {json: string, id: string} {
    const library = this.libraries.get(key);
    if (!library)
      throw new Error(`Chem | no monomer library ${key}: setLibrary first`);
    this.libraries.delete(key);
    this.libraries.set(key, library);
    return library;
  }

  private async mainThreadApi(): Promise<CruxBiopolymerApi> {
    const module = await loadCrux('biopolymer');
    if (this.compiledOf !== module) { // another instance (the tests' loader reset): its objects are not freed here
      this.compiled.clear();
      this.compiledOf = module;
    }
    return this.api = cruxBiopolymer(module);
  }

  // The main thread's module for a synchronous twin: loaded by an earlier ready() (or another async call), and not
  // after a trap there (Crux is off for the session, crux-module.ts noteCruxError: the instance may be broken, and the
  // async twins refuse it too)
  private syncApi(): CruxBiopolymerApi {
    if (!this.api) {
      throw new Error('Chem | the macromolecule assembler\'s synchronous calls need ready(key, json) first: ' +
        'the main thread\'s Crux module is not loaded');
    }
    if (cruxLoadFailed())
      throw new Error('Chem | the macromolecule assembler: the main thread\'s Crux module trapped in this session');
    return this.api;
  }

  private async mainThreadLibrary(id: string, json: string): Promise<CruxMonomerLibrary> {
    await this.mainThreadApi();
    return this.mainThreadLibrarySync(id, json);
  }

  // Compiled now where it is not (or was freed).
  private mainThreadLibrarySync(id: string, json: string): CruxMonomerLibrary {
    const crux = this.syncApi();
    let library = this.compiled.get(id);
    if (library)
      this.compiled.delete(id);
    else {
      library = onMainThread(() => new crux.ChemMonomerLibrary(json));
      while (this.compiled.size >= MAX_COMPILED)
        this.freeCompiled(this.compiled.keys().next().value!);
    }
    this.compiled.set(id, library);
    return library;
  }

  private freeCompiled(id: string): void {
    this.compiled.get(id)?.free();
    this.compiled.delete(id);
  }
}

/** The AssembleOneResult of a row Crux assembled; the map only when asked and the row has no error. */
function assembledOne(r: CruxAssembledRow, options: AssembleOptions): AssembleOneResult {
  const error = r.error ?? null;
  // the layout as Crux gives it: a value chem-meta's AsmLayout does not list yet is kept
  const result: AssembleOneResult = {molblock: error ? null : r.molblock ?? null, layout: r.layout as AsmLayout,
    error};
  if (options.map && !error && r.atomPosition && r.bondPosition)
    result.map = {atomPosition: int32s(r.atomPosition), bondPosition: int32s(r.bondPosition)};
  return result;
}

/** The MOL blocks markStereoAbsolute sends (non-empty strings) and the result with the empty cells in place. */
function stereoInput(molblocks: (string | null)[]): {result: (string | null)[], rows: number[], sent: string[]} {
  const result = Array.from(molblocks, (m): string | null => m === '' ? '' : null);
  const rows: number[] = [];
  const sent: string[] = [];
  Array.from(molblocks).forEach((m, i) => {
    if (typeof m === 'string' && m !== '') {
      rows.push(i);
      sent.push(m);
    }
  });
  return {result, rows, sent};
}

let assembler: CruxMacromoleculeAssembler | null = null;

/** The package's assembler (`Chem:getMacromoleculeAssembler`), one per session. */
export function getMacromoleculeAssembler(): CruxMacromoleculeAssembler {
  return assembler ??= new CruxMacromoleculeAssembler();
}
