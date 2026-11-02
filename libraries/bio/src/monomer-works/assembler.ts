import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {LRUCache} from 'lru-cache';

import {ChemTags} from '@datagrok-libraries/chem-meta/src/consts';
import {
  AsmAtomMap, AsmBatch, AsmLayout, AsmRole, AsmSequence, AssembleOneResult, AssembleResult, IMacromoleculeAssembler,
  unpackRow,
} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';

import {HelmType} from '../helm/types';
import {HelmTypes} from '../helm/consts';
import {IMonomerLibBase} from '../types/monomer-library';
import {ISeqHelper, ToAtomicLevelRes} from '../utils/seq-helper';
import {buildSeqBatch} from './assembly-input';
import {getMonomerLibJson, getMonomerLibKey} from './library-json';
import {MonomerMap} from './types';
import {getMolColName} from './utils';

/* Sequence to molecule through Chem's macromolecule assembler (Bio on Crux, PLAN §5.4). Every conversion, linear or
 * non-linear, column or single sequence, builds the parsed sequences (assembly-input.ts), registers the effective
 * library (library-json.ts) and calls the assembler. Column conversions return molblocks only; the monomer maps of
 * the hover are asked for per row ({@link MonomerMapProvider}). */

/** Column tag of a molecule column made from sequences: the layout it was assembled with, 'linear' or 'nonlinear'.
 * The hover asks for its maps with the same layout. */
export const SEQUENCE_SRC_LAYOUT_TAG = '.sequence-src-layout';

export type AsmLayoutName = 'linear' | 'nonlinear';

/** The layout a molecule column was assembled with: its {@link SEQUENCE_SRC_LAYOUT_TAG}, any value but 'linear'
 * read as 'nonlinear'; {@link defaultLayout} for columns without the tag. */
export function getSeqSrcLayout(molCol: DG.Column<string>, defaultLayout: AsmLayoutName): AsmLayoutName {
  const tag = molCol.getTag(SEQUENCE_SRC_LAYOUT_TAG);
  if (!tag) return defaultLayout;
  return tag === 'linear' ? 'linear' : 'nonlinear';
}

type AssemblerWindowType = Window & { $macromoleculeAssemblerPromise?: Promise<IMacromoleculeAssembler> };
declare const window: AssemblerWindowType;

/** The session's assembler (Chem loads Crux for it in both engine modes). A failed request is retried on the next
 * call. */
export async function getMacromoleculeAssembler(): Promise<IMacromoleculeAssembler> {
  let res = window.$macromoleculeAssemblerPromise;
  if (res === undefined) {
    const promise = (async (): Promise<IMacromoleculeAssembler> => {
      if (DG.Func.find({package: 'Chem', name: 'getMacromoleculeAssembler'}).length === 0) {
        throw new Error('Package "Chem" with the macromolecule assembler (Chem:getMacromoleculeAssembler) ' +
          'must be installed to convert sequences to molecules.');
      }
      const assembler = (await grok.functions.call('Chem:getMacromoleculeAssembler')) as IMacromoleculeAssembler;
      if (!assembler) throw new Error('Chem:getMacromoleculeAssembler returned no assembler.');
      return assembler;
    })();
    promise.catch(() => {
      if (window.$macromoleculeAssemblerPromise === promise) window.$macromoleculeAssemblerPromise = undefined;
    });
    res = window.$macromoleculeAssemblerPromise = promise;
  }
  return res;
}

/** Library keys registered per assembler object; the promise is dropped on failure so that it is sent again. */
const libRegistrations = new WeakMap<IMacromoleculeAssembler, Map<string, Promise<void>>>();

/** Registers the effective library with the assembler, once per library key.
 * @param {IMacromoleculeAssembler} assembler
 * @param {IMonomerLibBase} monomerLib the library to resolve with, an override library included
 * @return {Promise<string>} the library key to assemble with */
export async function setAssemblerLibrary(
  assembler: IMacromoleculeAssembler, monomerLib: IMonomerLibBase,
): Promise<string> {
  let byKey = libRegistrations.get(assembler);
  if (!byKey) libRegistrations.set(assembler, byKey = new Map<string, Promise<void>>());
  const key = getMonomerLibKey(monomerLib);
  let reg = byKey.get(key);
  if (!reg) {
    const lib = getMonomerLibJson(monomerLib);
    reg = assembler.setLibrary(lib.key, lib.json);
    byKey.set(key, reg);
    reg.catch(() => byKey!.delete(key));
  }
  await reg;
  return key;
}

/** Registers the effective library and readies it on the assembler's main-thread module ({@link
 * IMacromoleculeAssembler.ready}): after it, the synchronous twins can use the returned key.
 * @param {IMacromoleculeAssembler} assembler
 * @param {IMonomerLibBase} monomerLib the library to resolve with, an override library included
 * @return {Promise<string>} the library key for the synchronous calls */
export async function readyAssemblerLibrary(
  assembler: IMacromoleculeAssembler, monomerLib: IMonomerLibBase,
): Promise<string> {
  const lib = getMonomerLibJson(monomerLib);
  await assembler.ready(lib.key, lib.json);
  let byKey = libRegistrations.get(assembler);
  if (!byKey) libRegistrations.set(assembler, byKey = new Map<string, Promise<void>>());
  byKey.set(lib.key, Promise.resolve()); // ready() registers the library as setLibrary() does
  return lib.key;
}

/** Assembles every row of a macromolecule column (any notation) into a molecule column with the contract of the
 * old conversions: named `molfile(<seqCol>)` (an unused name in the column's table), Molecule / molblock, tags
 * `.sequence-src-col` and {@link SEQUENCE_SRC_LAYOUT_TAG}, V3000 molblocks, '' for empty and failed rows. Rows that
 * cannot be built go to warnings, one per row; more than 5 % of them show a warning when {@link notify}.
 * The column is not added to the table.
 * @param {DG.Column<string>} seqCol a Macromolecule column
 * @param {IMonomerLibBase} monomerLib the library to resolve with
 * @param {ISeqHelper} seqHelper
 * @param {AsmLayoutName} layout 'linear' falls back per row to the non-linear layout
 * @param {boolean} notify the warning balloon; off for one-sequence callers (the panels show their own message)
 * @return {Promise<ToAtomicLevelRes>} */
export async function assembleSeqColumn(
  seqCol: DG.Column<string>, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, layout: AsmLayoutName,
  notify: boolean = true,
): Promise<ToAtomicLevelRes> {
  const notMacromolecule = checkMacromolecule(seqCol);
  if (notMacromolecule) return notMacromolecule;
  const batch = buildSeqBatch(seqHelper.getSeqHandler(seqCol), null, monomerLib);
  const assembler = await getMacromoleculeAssembler();
  const libKey = await setAssemblerLibrary(assembler, monomerLib);
  return columnOfResult(seqCol, batch, await assembler.assemble(batch, libKey, {layout}), layout, notify);
}

/** {@link assembleSeqColumn} on the assembler's main-thread module, synchronously (the legacy synchronous APIs):
 * {@link libKey} must have been readied ({@link readyAssemblerLibrary}).
 * @return {ToAtomicLevelRes} */
export function assembleSeqColumnSync(
  seqCol: DG.Column<string>, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper, layout: AsmLayoutName,
  assembler: IMacromoleculeAssembler, libKey: string,
): ToAtomicLevelRes {
  const notMacromolecule = checkMacromolecule(seqCol);
  if (notMacromolecule) return notMacromolecule;
  const batch = buildSeqBatch(seqHelper.getSeqHandler(seqCol), null, monomerLib);
  return columnOfResult(seqCol, batch, assembler.assembleSync(batch, libKey, {layout}), layout);
}

function checkMacromolecule(seqCol: DG.Column<string>): ToAtomicLevelRes | null {
  if (seqCol.semType === DG.SEMTYPE.MACROMOLECULE) return null;
  const msg: string = `Only the ${DG.SEMTYPE.MACROMOLECULE} columns can be converted to atomic level, ` +
    `the chosen column has semType '${seqCol.semType}'`;
  return {molCol: null, warnings: [msg]};
}

/** The result column of {@link assembleSeqColumn}: per-row errors to warnings, the toast when {@link notify}, the
 * tags. */
function columnOfResult(
  seqCol: DG.Column<string>, batch: AsmBatch, res: AssembleResult, layout: AsmLayoutName, notify: boolean = true,
): ToAtomicLevelRes {
  const rowCount = seqCol.length;
  const warnings: string[] = [];
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const err = res.errors[rowIdx] ?? batch.rowError[rowIdx];
    if (err) warnings.push(`Cannot get molfile of row #${rowIdx}: ${err}.`);
  }
  if (notify && warnings.length > 0.05 * rowCount)
    grok.shell.warning(`Molfile conversion resulted in ${warnings.length} errors`);

  const molColName = getMolColName(seqCol.dataFrame, seqCol.name);
  const molCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, molColName, rowCount)
    .init((rowIdx) => res.molblocks[rowIdx] ?? '') as DG.Column<string>;
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;
  molCol.setTag(ChemTags.SEQUENCE_SRC_COL, seqCol.name);
  molCol.setTag(SEQUENCE_SRC_LAYOUT_TAG, layout);
  return {molCol, warnings};
}

function biotypeOfRole(role: AsmRole): HelmType {
  switch (role) {
  case AsmRole.AA: return HelmTypes.AA as HelmType;
  case AsmRole.CHEM: return HelmTypes.CHEM as HelmType;
  default: return HelmTypes.NUCLEOTIDE as HelmType;
  }
}

/** Bio's monomer map of one assembled row: per position (the map key), the atoms and bonds it owns, 0-based in
 * the row's molblock. Several monomers share a position (the copies of an expanded repeat; a FASTA/separator
 * nucleotide's phosphate, sugar and base): their atoms and bonds are merged, the symbol and biotype come from the
 * base if there is one, else from the first monomer. RNA monomers get HELM_NUCLEOTIDE as in the old maps (the
 * biotype getMolHighlight looks the colour up with); an inline-SMILES monomer's symbol is its SMILES.
 * @param {AsmSequence} seq the assembled sequence
 * @param {AsmAtomMap} map the assembler's provenance of the row
 * @return {MonomerMap} */
export function asmMapToMonomerMap(seq: AsmSequence, map: AsmAtomMap): MonomerMap {
  const res = new MonomerMap();
  const roleByPos = new Map<number, AsmRole>();
  for (const m of seq.monomers) {
    if (m.gap) continue;
    const prevRole = roleByPos.get(m.position);
    if (prevRole !== undefined && (prevRole === AsmRole.BASE || m.role !== AsmRole.BASE)) continue;
    roleByPos.set(m.position, m.role);
    res.set(m.position, {biotype: biotypeOfRole(m.role), symbol: m.symbol || m.inlineSmiles || '',
      atoms: [], bonds: []});
  }
  for (let atomIdx = 0; atomIdx < map.atomPosition.length; ++atomIdx) {
    const pos = map.atomPosition[atomIdx];
    if (pos >= 0) res.get(pos)?.atoms.push(atomIdx);
  }
  for (let bondIdx = 0; bondIdx < map.bondPosition.length; ++bondIdx) {
    const pos = map.bondPosition[bondIdx];
    if (pos >= 0) res.get(pos)?.bonds.push(bondIdx);
  }
  return res;
}

/** The parsed sequence of one row, or null for an empty or rejected row. */
export function getRowAsmSequence(
  seqCol: DG.Column<string>, rowIdx: number, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper,
): { seq: AsmSequence | null, error: string | null } {
  const batch = buildSeqBatch(seqHelper.getSeqHandler(seqCol), [rowIdx], monomerLib);
  const error = batch.rowError[0];
  if (error) return {seq: null, error};
  const seq = unpackRow(batch, 0);
  return {seq: seq.monomers.length > 0 ? seq : null, error: null};
}

/** One row assembled on the main-thread module with its monomer map ('' and an empty map for an empty or failed
 * row). Used by the converter's wrapper and the hover.
 * @return {Promise<AssembledRow>} */
export async function assembleRowWithMap(
  seqCol: DG.Column<string>, rowIdx: number, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper,
  layout: AsmLayoutName,
): Promise<AssembledRow> {
  const {seq, error} = getRowAsmSequence(seqCol, rowIdx, monomerLib, seqHelper);
  if (!seq) return {molfile: '', monomers: new MonomerMap(), layout: null, error};
  const assembler = await getMacromoleculeAssembler();
  const libKey = await setAssemblerLibrary(assembler, monomerLib);
  return rowOfResult(seq, await assembler.assembleOne(seq, libKey, {layout, map: true}));
}

/** {@link assembleRowWithMap} on the assembler's main-thread module, synchronously (the legacy synchronous APIs):
 * {@link libKey} must have been readied ({@link readyAssemblerLibrary}). */
export function assembleRowWithMapSync(
  seqCol: DG.Column<string>, rowIdx: number, monomerLib: IMonomerLibBase, seqHelper: ISeqHelper,
  layout: AsmLayoutName, assembler: IMacromoleculeAssembler, libKey: string,
): AssembledRow {
  const {seq, error} = getRowAsmSequence(seqCol, rowIdx, monomerLib, seqHelper);
  if (!seq) return {molfile: '', monomers: new MonomerMap(), layout: null, error};
  return rowOfResult(seq, assembler.assembleOneSync(seq, libKey, {layout, map: true}));
}

export type AssembledRow = { molfile: string, monomers: MonomerMap, layout: AsmLayout | null, error: string | null };

function rowOfResult(seq: AsmSequence, res: AssembleOneResult): AssembledRow {
  if (res.error || !res.molblock)
    return {molfile: '', monomers: new MonomerMap(), layout: res.layout, error: res.error ?? 'no molecule'};
  const monomers = res.map ? asmMapToMonomerMap(seq, res.map) : new MonomerMap();
  return {molfile: res.molblock, monomers, layout: res.layout, error: null};
}

const MAP_CACHE_MAX = 100;
/** Requests waiting beyond this count are dropped, the oldest first (rows scrolled away while their maps waited) */
const MAP_QUEUE_MAX = 200;

type MapRequest = { key: string, rowIdx: number, seq: string, layout: AsmLayoutName };
type MapReady = () => void;

/** The monomer maps of a molecule column made from a sequence column, from the assembler on demand, one row at a
 * time (D6). The renderers ask synchronously: {@link get} returns the map when it is cached and otherwise requests
 * it and returns undefined; the callback runs when the map arrives. Maps are cached by layout and sequence string
 * (the last {@link MAP_CACHE_MAX}); the requests run one after another, the hovered row first. */
export class MonomerMapProvider {
  /** false: the row has no map (empty, rejected or failed) */
  private readonly cache = new LRUCache<string, MonomerMap | false>({max: MAP_CACHE_MAX});
  private readonly waiting = new Map<string, MapReady[]>();
  private readonly queue: MapRequest[] = [];
  private running: boolean = false;

  constructor(
    public readonly seqCol: DG.Column<string>,
    private readonly monomerLib: IMonomerLibBase,
    private readonly seqHelper: ISeqHelper,
    private readonly getLayout: (tableRowIdx: number) => AsmLayoutName,
  ) {}

  /** @param {number} tableRowIdx
   * @param {Function} onReady called when a requested map arrives (not when it is cached already)
   * @param {boolean} urgent put the request first (the hovered row)
   * @return {MonomerMap | null | undefined} the map; null when the row has none; undefined while it is requested */
  get(tableRowIdx: number, onReady?: MapReady, urgent: boolean = false): MonomerMap | null | undefined {
    const seq = this.seqCol.get(tableRowIdx);
    if (!seq) return null;
    const layout = this.getLayout(tableRowIdx);
    const key = `${layout}\u0001${seq}`;
    const cached = this.cache.get(key);
    if (cached !== undefined) return cached || null;

    let callbacks = this.waiting.get(key);
    if (!callbacks) {
      this.waiting.set(key, callbacks = []);
      const req: MapRequest = {key, rowIdx: tableRowIdx, seq, layout};
      if (urgent) this.queue.unshift(req); else this.queue.push(req);
      while (this.queue.length > MAP_QUEUE_MAX) {
        const dropIdx = this.queue.findIndex((r) => r.key !== key);
        const [dropped] = this.queue.splice(dropIdx < 0 ? 0 : dropIdx, 1);
        this.waiting.delete(dropped.key);
      }
      this.run();
    } else if (urgent) {
      const reqIdx = this.queue.findIndex((r) => r.key === key);
      if (reqIdx > 0) this.queue.unshift(...this.queue.splice(reqIdx, 1));
    }
    if (onReady && !callbacks.includes(onReady)) callbacks.push(onReady);
    return undefined;
  }

  private async run(): Promise<void> {
    if (this.running) return;
    this.running = true;
    try {
      while (this.queue.length > 0) {
        const req = this.queue.shift()!;
        let map: MonomerMap | false = false;
        let current = true;
        try {
          // The row may have changed while its request waited: no map is cached for a stale request
          current = this.seqCol.get(req.rowIdx) === req.seq;
          if (current) {
            const res = await assembleRowWithMap(this.seqCol, req.rowIdx, this.monomerLib, this.seqHelper,
              req.layout);
            map = res.error ? false : res.monomers;
          }
        } catch (err: any) {
          console.error(err);
        }
        if (current) this.cache.set(req.key, map);
        const callbacks = this.waiting.get(req.key) ?? [];
        this.waiting.delete(req.key);
        for (const cb of callbacks) {
          try {
            cb();
          } catch (err: any) {
            console.error(err);
          }
        }
      }
    } finally {
      this.running = false;
    }
  }
}
