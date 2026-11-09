/**
 * The Crux worker's macromolecule assembly (Bio-on-Crux PLAN §5.3; crux.worker.ts runs it, CruxService
 * .assembleMacromolecules sends it): a monomer library is compiled once per worker (`biopolymer::library`, kept by
 * key, the last {@link BIOPOLYMER_LIBRARIES} used, the least recently used freed first, as jobs/batch.ts keeps its
 * SMARTS sets), then each chunk of a packed batch (chem-meta's AsmBatch) is assembled in one call
 * (`biopolymer::assemble`). A batch naming a library the worker does not hold fails with
 * {@link BIOPOLYMER_LIBRARY_MISSING}: the service sends the library and asks once more.
 */
import type {AsmBatch} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import {BIOPOLYMER_LIBRARIES, BIOPOLYMER_LIBRARY_MISSING, CruxAssembledBatch, CruxMonomerLibrary, WorkerJob,
  cruxBiopolymer, throwIfTrap} from '../crux-api';

/** The libraries by key, the most recently used last, each with the number of its monomers Crux does not read. */
const libraries = new Map<string, {library: CruxMonomerLibrary, errors: number}>();
/** Tests: the libraries this worker has compiled. */
let compiledCount = 0;

/** The number of monomers `library` cannot read (its `errors()` JSON, a list; read once: it compiles every
 * template). */
function readErrorCount(library: CruxMonomerLibrary): number {
  try {
    const text = library.errors();
    const errors = text ? JSON.parse(text) : null;
    return Array.isArray(errors) ? errors.length : 0;
  } catch (e) {
    throwIfTrap(e);
    return 0; // only counted (a log line): no reason to refuse the library
  }
}

/** The library held under `key` (the most recently used now), or the error the service knows. */
function libraryOf(key: string): {library: CruxMonomerLibrary, errors: number} {
  const held = libraries.get(key);
  if (!held)
    throw new Error(`${BIOPOLYMER_LIBRARY_MISSING} ${key}`);
  libraries.delete(key);
  libraries.set(key, held);
  return held;
}

/** Compiles `json` (Pistoia monomer library JSON) under `key`, unless it is held already: the number of monomers it
 * does not read. */
function loadLibrary(crux: unknown, key: string, json: string): {errors: number} {
  if (libraries.has(key))
    return {errors: libraryOf(key).errors};
  const {ChemMonomerLibrary} = cruxBiopolymer(crux);
  let library: CruxMonomerLibrary;
  try {
    library = new ChemMonomerLibrary(json);
  } catch (e: any) {
    throwIfTrap(e);
    throw new Error(`biopolymer: the monomer library ${key} does not read: ${e?.message ?? e}`);
  }
  compiledCount++;
  while (libraries.size >= BIOPOLYMER_LIBRARIES) {
    const [oldest, entry] = libraries.entries().next().value!;
    entry.library.free();
    libraries.delete(oldest);
  }
  const errors = readErrorCount(library);
  libraries.set(key, {library, errors});
  return {errors};
}

/** `array` itself where it owns its buffer (transferable alone), else a copy that does (a view into the wasm memory
 * must not be transferred). */
function ownBuffer(array: Uint8Array): Uint8Array {
  return array.buffer instanceof ArrayBuffer && array.byteOffset === 0 &&
    array.byteLength === array.buffer.byteLength ? array : array.slice();
}

/**
 * The rows of `batch` assembled with library `key` (`options`: AssembleOptions as JSON): a MOL block (null where the
 * row has an error), a layout and an error per row. The layouts' buffer is transferred.
 */
function assemble(transfer: Transferable[], batch: AsmBatch, key: string, options: string): CruxAssembledBatch {
  const {library} = libraryOf(key);
  const n = batch.rowStart.length - 1;
  const r = library.assemble(JSON.stringify(batch.symbols), batch.rowStart, batch.monSymbol, batch.monRole,
    batch.monPosition, batch.chainStart, batch.rowChainStart, batch.conn, batch.rowConnStart,
    Array.from(batch.rowError, (e) => e ?? null), options);
  const errors = Array.from({length: n}, (_, i) => r.errors[i] ?? null);
  const molblocks = Array.from({length: n}, (_, i) => errors[i] ? null : r.molblocks[i] ?? null);
  const layouts = ownBuffer(r.layouts);
  transfer.push(layouts.buffer);
  return {molblocks, layouts, errors};
}

/** The jobs by name (the Crux worker's `op`). */
export const BIOPOLYMER_JOBS: {[op: string]: WorkerJob} = {
  'biopolymer::library': ({crux}, key: string, json: string) => loadLibrary(crux, key, json),
  'biopolymer::assemble': ({transfer}, batch: AsmBatch, key: string, options: string) =>
    assemble(transfer, batch, key, options),
  /** A column's MOL blocks with their stereocentres in one ABS group (engine/biopolymer.ts, a column's call). */
  'biopolymer::markStereoAbsolute': ({crux}, molblocks: (string | null)[]) =>
    cruxBiopolymer(crux).chemMarkStereoAbsolute(molblocks),
  /** Tests: the keys this worker holds (the least recently used first) and how many libraries it has compiled. */
  'biopolymer::libraries::test': () => ({keys: [...libraries.keys()], compiled: compiledCount}),
  /** Tests: the worker forgets library `key` (its client still counts it as sent). */
  'biopolymer::forget::test': (_, key: string) => {
    libraries.get(key)?.library.free();
    libraries.delete(key);
  },
};
