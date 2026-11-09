/** Interfaces of the macromolecule assembler (Chem:getMacromoleculeAssembler) and its input.
 *
 * Datagrok parses every notation and sends a notation-free parsed sequence ({@link AsmSequence}): monomers with
 * exact library keys, chains, explicit connections. The assembler builds the molecule from the monomer library's
 * templates. Implicit bonds are not sent; the assembler derives them:
 * - within a chain, consecutive non-gap backbone monomers: R2 of the left one to R1 of the right one;
 * - RNA: sugar R3 to base R1 (the base follows its sugar), sugar R2 to linker R1, linker R2 to the next sugar R1;
 * - a CHEM chain is one monomer with no implicit bonds.
 *
 * {@link AsmSequence} is the readable form (tests, single-row calls); {@link AsmBatch} is the packed form sent to
 * workers ({@link packBatch} / {@link unpackRow}).
 */

export type AsmPolymerType = 'PEPTIDE' | 'RNA' | 'CHEM' | 'BLOB';

/** HelmType equivalents: AA, SUGAR, BASE, LINKER, CHEM (BLOB rows never reach the assembler). */
export const enum AsmRole { AA = 0, SUGAR = 1, BASE = 2, LINKER = 3, CHEM = 4 }

export interface AsmMonomer {
  /** Exact library symbol (after Datagrok's resolution), or '' for an inline-SMILES monomer or a gap. */
  symbol: string;
  /** Exact library polymer type (with the symbol, the library key); for a gap, the chain's polymer type. */
  polymerType: AsmPolymerType;
  role: AsmRole;
  /** Inline-SMILES monomer: R-groups as [*:n] or CXSMILES _Rn labels. */
  inlineSmiles?: string;
  /** Gap ('*' in HELM, '-' or '' elsewhere): skipped, the backbone continues across it. */
  gap?: boolean;
  /** Map key: the position in the original sequence (Datagrok's flat index; HELM RNA counts sugar, base and
   * linker). Equal to the monomer's index except for expanded repeats, whose copies keep the source token's
   * position, and for FASTA/separator nucleotides, whose sugar, base and phosphate all carry the base's index. */
  position: number;
}

/** Datagrok's ISeqConnection plus kind; indices 0-based, gaps counted. */
export interface AsmConnection {
  /** Chain index in source order */
  seqIndex1: number;
  /** Monomer index within the chain */
  monomerIndex1: number;
  /** R-group number, any integer (R10+ allowed) */
  rGroup1: number;
  seqIndex2: number;
  monomerIndex2: number;
  rGroup2: number;
  /** hbond: no bond in the molecule; default covalent */
  kind?: 'covalent' | 'hbond';
}

export interface AsmGraphInfo {
  /** Explicit connections only (HELM section 1, BILN marks); backbone and branch bonds are implicit */
  connections: AsmConnection[];
  /** Chain starts in `monomers` ([0] === 0 when there are monomers) */
  disjointSeqStarts: number[];
  /** 'PEPTIDE1', 'RNA2', ...; diagnostics only, not packed */
  chainIds?: string[];
}

export interface AsmSequence {
  /** Flat, chain after chain */
  monomers: AsmMonomer[];
  graphInfo: AsmGraphInfo;
}

/** An entry of the batch's symbol table. */
export interface AsmSymbol {
  polymerType: AsmPolymerType;
  symbol: string;
  inlineSmiles?: string;
}

/** Number of Int32 values per connection in {@link AsmBatch.conn}: s1, m1, r1, s2, m2, r2, kind. */
export const ASM_CONN_STRIDE = 7;

/** Connection kind codes in {@link AsmBatch.conn}. */
export const enum AsmConnKind { COVALENT = 0, HBOND = 1 }

/** The packed batch, sent per chunk; its typed arrays are transferable ({@link asmBatchTransferables}). */
export interface AsmBatch {
  /** Distinct per batch */
  symbols: AsmSymbol[];
  /** n+1 offsets into the monomer arrays */
  rowStart: Int32Array;
  /** Index into symbols; -1 = gap */
  monSymbol: Int32Array;
  /** AsmRole */
  monRole: Uint8Array;
  /** Map key */
  monPosition: Int32Array;
  /** Flattened disjointSeqStarts (row-local) */
  chainStart: Int32Array;
  /** n+1 offsets into chainStart */
  rowChainStart: Int32Array;
  /** {@link ASM_CONN_STRIDE} values per connection: s1, m1, r1, s2, m2, r2, kind (0 covalent, 1 hbond) */
  conn: Int32Array;
  /** n+1 offsets into conn, in connections (multiply by {@link ASM_CONN_STRIDE} for the Int32 offset) */
  rowConnStart: Int32Array;
  /** Rows Datagrok already rejected (unknown monomer, ambiguity, BLOB, range...); they carry no monomers and the
   * assembler skips them */
  rowError: (string | null)[];
}

export interface AssembleOptions {
  /** linear falls back per row to the non-linear layout */
  layout: 'linear' | 'nonlinear';
  /** Default 'depictor'. 'auto': CoordGen for big rows (clean drawings of long chains, about 8x slower), the depictor
   * for the rest */
  nonlinearEngine?: 'depictor' | 'coordgen' | 'auto';
  /** Single-row calls only */
  map?: boolean;
}

/** How a row was laid out: 0 linear; 1 non-linear (depictor); 2 linear asked, fallen back to non-linear (depictor);
 * 3 non-linear drawn by CoordGen (asked, or the assembler's size switch for big rows); 4 linear asked, fallen back
 * and drawn by CoordGen; 255 error. Any value but LINEAR and ERROR is a non-linear layout (more may be added). */
export const enum AsmLayout {
  LINEAR = 0, NONLINEAR = 1, NONLINEAR_FALLBACK = 2, NONLINEAR_COORDGEN = 3, NONLINEAR_FALLBACK_COORDGEN = 4,
  ERROR = 255,
}

export interface AssembleResult {
  /** V3000; null where errors[i] */
  molblocks: (string | null)[];
  /** AsmLayout per row */
  layouts: Uint8Array;
  /** Per-row message ("unknown monomer PEPTIDE:Xyz", "R-group R3 of ... not found", ...) */
  errors: (string | null)[];
}

/** One row's atom and bond provenance. */
export interface AsmAtomMap {
  /** Per atom of the molblock: the monomer's position (map key), -1 none */
  atomPosition: Int32Array;
  /** Per bond: owning position, -1 none */
  bondPosition: Int32Array;
}

export interface AssembleOneResult {
  molblock: string | null;
  layout: AsmLayout;
  error: string | null;
  map?: AsmAtomMap;
}

export interface IMacromoleculeAssembler {
  /** Registers a library (Pistoia JSON array text) under a key; idempotent; the key is the cache key in every
   * worker. */
  setLibrary(key: string, json: string): Promise<void>;

  assemble(batch: AsmBatch, libraryKey: string, options: AssembleOptions): Promise<AssembleResult>;

  /** One sequence on the main-thread module (panel, hover): molblock and, with options.map, the map. */
  assembleOne(seq: AsmSequence, libraryKey: string, options: AssembleOptions): Promise<AssembleOneResult>;

  /** Every specified stereocentre into one ABS group, V3000 out: replaces Chem:convertToV3KViaOCL for Bio's
   * callers. */
  markStereoAbsolute(molblocks: (string | null)[]): Promise<(string | null)[]>;

  /** V3000 to V2000 (Bio's 3D panel). */
  toV2000(molblock: string): Promise<string>;

  // -- Synchronous twins, for the legacy synchronous APIs (Bio's IHelmToMolfileConverter, helmToAtomicLevelSingle):
  // all on the main-thread module, so they work only after ready() and throw a clear error otherwise. New code uses
  // the async methods above.

  /** Registers the library as {@link setLibrary} does, loads the main-thread module and compiles the library there:
   * after it, the synchronous twins can use the key. */
  ready(libraryKey: string, json: string): Promise<void>;

  /** {@link assembleOne} on the main-thread module, synchronously; the key must have been readied. */
  assembleOneSync(seq: AsmSequence, libraryKey: string, options: AssembleOptions): AssembleOneResult;

  /** {@link assemble} on the main-thread module (not the workers), synchronously; the key must have been readied. */
  assembleSync(batch: AsmBatch, libraryKey: string, options: AssembleOptions): AssembleResult;

  /** {@link markStereoAbsolute} on the main-thread module, synchronously; after any ready(). */
  markStereoAbsoluteSync(molblocks: (string | null)[]): (string | null)[];

  /** {@link toV2000} on the main-thread module, synchronously; after any ready(). */
  toV2000Sync(molblock: string): string;
}

const POLYMER_TYPE_OF_ROLE: AsmPolymerType[] = ['PEPTIDE', 'RNA', 'RNA', 'RNA', 'CHEM'];

/** The polymer type a gap takes when unpacked: gaps carry no symbol, so it follows the role. */
export function asmPolymerTypeOfRole(role: AsmRole): AsmPolymerType {
  return POLYMER_TYPE_OF_ROLE[role] ?? 'PEPTIDE';
}

/** Checks the invariants {@link packBatch} relies on.
 * @return {string | null} the first problem found, or null */
export function validateAsmSequence(seq: AsmSequence): string | null {
  const n = seq.monomers.length;
  const starts = seq.graphInfo.disjointSeqStarts;
  if (n > 0 && (starts.length === 0 || starts[0] !== 0))
    return 'disjointSeqStarts must start with 0';
  for (let i = 0; i < starts.length; ++i) {
    if (!Number.isInteger(starts[i]) || starts[i] < 0 || starts[i] >= Math.max(n, 1) ||
      (i > 0 && starts[i] <= starts[i - 1]))
      return `invalid chain start ${starts[i]} at ${i}`;
  }
  for (let i = 0; i < n; ++i) {
    const m = seq.monomers[i];
    if (!m.gap && !m.symbol && !m.inlineSmiles)
      return `monomer ${i} has neither a symbol nor inline SMILES`;
    if (!Number.isInteger(m.position))
      return `monomer ${i} has no position`;
  }
  const chainLength = (s: number): number => (s + 1 < starts.length ? starts[s + 1] : n) - starts[s];
  for (const c of seq.graphInfo.connections) {
    for (const [s, m, r] of [[c.seqIndex1, c.monomerIndex1, c.rGroup1], [c.seqIndex2, c.monomerIndex2, c.rGroup2]]) {
      if (!Number.isInteger(s) || s < 0 || s >= starts.length)
        return `connection chain index ${s} out of range`;
      if (!Number.isInteger(m) || m < 0 || m >= chainLength(s))
        return `connection monomer index ${m} out of range of chain ${s}`;
      if (!Number.isInteger(r) || r < 0)
        return `connection R-group ${r} is invalid`;
    }
  }
  return null;
}

function symbolKey(polymerType: string, symbol: string, inlineSmiles: string | undefined): string {
  return `${polymerType}\u0001${symbol}\u0001${inlineSmiles ?? ''}`;
}

/** Packs sequences into a batch. A row is an error row when {@link rowErrors}[i] is set or {@link seqs}[i] is null;
 * error rows carry no monomers. Throws on a sequence that breaks {@link validateAsmSequence}.
 * @param {(AsmSequence | null)[]} seqs one entry per row
 * @param {(string | null)[]} rowErrors optional per-row errors
 * @return {AsmBatch} */
export function packBatch(seqs: (AsmSequence | null)[], rowErrors?: (string | null)[]): AsmBatch {
  const rowCount = seqs.length;
  const rowError: (string | null)[] = new Array(rowCount);
  let monCount = 0;
  let chainCount = 0;
  let connCount = 0;
  for (let i = 0; i < rowCount; ++i) {
    const seq = seqs[i];
    const err = rowErrors?.[i] ?? (seq ? null : 'no sequence');
    rowError[i] = err;
    if (err || !seq) continue;
    const invalid = validateAsmSequence(seq);
    if (invalid)
      throw new Error(`Row ${i}: ${invalid}`);
    monCount += seq.monomers.length;
    chainCount += seq.graphInfo.disjointSeqStarts.length;
    connCount += seq.graphInfo.connections.length;
  }

  const symbols: AsmSymbol[] = [];
  const symbolIdx = new Map<string, number>();
  const rowStart = new Int32Array(rowCount + 1);
  const monSymbol = new Int32Array(monCount);
  const monRole = new Uint8Array(monCount);
  const monPosition = new Int32Array(monCount);
  const rowChainStart = new Int32Array(rowCount + 1);
  const chainStart = new Int32Array(chainCount);
  const rowConnStart = new Int32Array(rowCount + 1);
  const conn = new Int32Array(connCount * ASM_CONN_STRIDE);

  let mi = 0;
  let ci = 0;
  let ki = 0;
  for (let i = 0; i < rowCount; ++i) {
    const seq = seqs[i];
    if (!rowError[i] && seq) {
      for (const m of seq.monomers) {
        if (m.gap)
          monSymbol[mi] = -1;
        else {
          const key = symbolKey(m.polymerType, m.symbol, m.inlineSmiles);
          let si = symbolIdx.get(key);
          if (si === undefined) {
            si = symbols.length;
            symbolIdx.set(key, si);
            const s: AsmSymbol = {polymerType: m.polymerType, symbol: m.symbol};
            if (m.inlineSmiles) s.inlineSmiles = m.inlineSmiles;
            symbols.push(s);
          }
          monSymbol[mi] = si;
        }
        monRole[mi] = m.role;
        monPosition[mi] = m.position;
        ++mi;
      }
      for (const s of seq.graphInfo.disjointSeqStarts)
        chainStart[ci++] = s;
      for (const c of seq.graphInfo.connections) {
        const o = ki * ASM_CONN_STRIDE;
        conn[o] = c.seqIndex1;
        conn[o + 1] = c.monomerIndex1;
        conn[o + 2] = c.rGroup1;
        conn[o + 3] = c.seqIndex2;
        conn[o + 4] = c.monomerIndex2;
        conn[o + 5] = c.rGroup2;
        conn[o + 6] = c.kind === 'hbond' ? AsmConnKind.HBOND : AsmConnKind.COVALENT;
        ++ki;
      }
    }
    rowStart[i + 1] = mi;
    rowChainStart[i + 1] = ci;
    rowConnStart[i + 1] = ki;
  }
  return {symbols, rowStart, monSymbol, monRole, monPosition, chainStart, rowChainStart, conn, rowConnStart,
    rowError};
}

/** Unpacks one row of a batch. An error row unpacks to an empty sequence; check {@link AsmBatch.rowError}.
 * Gaps get the polymer type of their role ({@link asmPolymerTypeOfRole}); chainIds are not packed.
 * @param {AsmBatch} batch
 * @param {number} i row index
 * @return {AsmSequence} */
export function unpackRow(batch: AsmBatch, i: number): AsmSequence {
  if (i < 0 || i >= batch.rowError.length)
    throw new Error(`Row ${i} is out of the batch's ${batch.rowError.length} rows`);
  const monomers: AsmMonomer[] = [];
  for (let mi = batch.rowStart[i]; mi < batch.rowStart[i + 1]; ++mi) {
    const role = batch.monRole[mi] as AsmRole;
    const si = batch.monSymbol[mi];
    if (si < 0) {
      const polymerType = asmPolymerTypeOfRole(role);
      monomers.push({symbol: '', polymerType, role, gap: true, position: batch.monPosition[mi]});
    } else {
      const s = batch.symbols[si];
      const m: AsmMonomer = {symbol: s.symbol, polymerType: s.polymerType, role, position: batch.monPosition[mi]};
      if (s.inlineSmiles) m.inlineSmiles = s.inlineSmiles;
      monomers.push(m);
    }
  }
  const disjointSeqStarts = Array.from(batch.chainStart.subarray(batch.rowChainStart[i], batch.rowChainStart[i + 1]));
  const connections: AsmConnection[] = [];
  for (let ki = batch.rowConnStart[i]; ki < batch.rowConnStart[i + 1]; ++ki) {
    const o = ki * ASM_CONN_STRIDE;
    const c: AsmConnection = {
      seqIndex1: batch.conn[o], monomerIndex1: batch.conn[o + 1], rGroup1: batch.conn[o + 2],
      seqIndex2: batch.conn[o + 3], monomerIndex2: batch.conn[o + 4], rGroup2: batch.conn[o + 5],
    };
    if (batch.conn[o + 6] === AsmConnKind.HBOND) c.kind = 'hbond';
    connections.push(c);
  }
  return {monomers, graphInfo: {connections, disjointSeqStarts}};
}

/** The batch's buffers, for postMessage's transfer list. */
export function asmBatchTransferables(batch: AsmBatch): ArrayBuffer[] {
  return [batch.rowStart, batch.monSymbol, batch.monRole, batch.monPosition, batch.chainStart, batch.rowChainStart,
    batch.conn, batch.rowConnStart].map((a) => a.buffer as ArrayBuffer);
}
