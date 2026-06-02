/*
 * Local (user-supplied) PDB registry.
 *
 * Supports dragging a `.pdb` file from the OS onto the Step 1 PDB-IDs textarea.
 * RCSB-coded files (e.g. `1XKK.pdb`) are added as normal IDs and flow through
 * the usual RCSB fetch + QC path. A *custom* file (e.g. `my_kinase.pdb`) has no
 * RCSB entry, so its raw text is registered here under a sanitized name and the
 * pipeline is routed around RCSB:
 *   - `fetchPdbBlock` (pdb-utils) returns the stored text instead of downloading,
 *   - `ensurePdbQc` (orchestrator) synthesizes a QC row from `buildLocalQcRows`
 *     instead of calling `enrichPdbList`,
 *   - `parsePdbIdsDetailed` (orchestrator) accepts the registered name even
 *     though it doesn't match the 4-char RCSB pattern.
 *
 * The drug ligand is auto-detected as the largest non-solvent HETATM residue
 * (see `detectLigand`). SMILES/MW are left empty for local structures (we can't
 * reliably perceive bond orders from coordinates), so the MW filter is skipped;
 * the pocket / feature stages locate the ligand by comp-id + chain + coords.
 */

/** RCSB 4-char code: a digit followed by three alphanumerics. */
const RCSB_ID_RE = /^[0-9][A-Z0-9]{3}$/;

/** Solvent / ion / common-buffer HETATM residues that are never the ligand.
 *  Intentionally local to this module so it has no dependency on pdb-utils
 *  (which depends on this module via fetchPdbBlock — avoids an import cycle). */
const NON_LIGAND_HET = new Set([
  'HOH', 'DOD', 'SO4', 'PO4', 'CL', 'NA', 'K', 'CA', 'MG', 'ZN', 'FE', 'MN',
  'NI', 'CU', 'EDO', 'GOL', 'PEG', 'PG4', 'PGE', 'DMS', 'ACT', 'TRS', 'MES',
  'IMD', 'NAG', 'BMA', 'MAN', 'FUC', 'GAL',
]);

export interface LocalPdb {
  /** Canonical id used in the textarea + pipeline (uppercase, [A-Z0-9_]). */
  id: string;
  /** Original dropped file name (for messages). */
  fileName: string;
  /** Raw PDB text. */
  pdbText: string;
  /** Parsed from `REMARK   2 RESOLUTION`, else null. */
  resolution: number | null;
  /** Auto-detected drug-ligand 3-letter residue name ('' if none found). */
  ligandCompId: string;
  /** Chain the detected ligand sits on ('' if none). */
  ligandChain: string;
  /** Cheap content hash, for cache fingerprinting. */
  contentHash: string;
}

const store = new Map<string, LocalPdb>();

export function isValidRcsbId(token: string): boolean {
  return RCSB_ID_RE.test(token.toUpperCase());
}

export function isLocalId(id: string): boolean {
  return store.has(id.toUpperCase());
}

export function getLocalPdb(id: string): LocalPdb | undefined {
  return store.get(id.toUpperCase());
}

export function localIdsIn(ids: string[]): string[] {
  return ids.filter((id) => store.has(id.toUpperCase()));
}

/** Sanitize a file name into a unique, textarea-safe id (uppercase, no
 *  whitespace/comma/semicolon so it survives the textarea split). */
function makeLocalId(fileName: string): string {
  const stem = fileName.replace(/\.[^.]*$/, '');           // drop extension
  let base = stem.toUpperCase().replace(/[^A-Z0-9]+/g, '_')
    .replace(/^_+|_+$/g, '');
  if (!base) base = 'LOCAL';
  if (!/[A-Z]/.test(base)) base = 'L_' + base;             // keep it un-RCSB-like
  let id = base;
  let n = 2;
  while (store.has(id)) id = `${base}_${n++}`;             // de-dupe
  return id;
}

function hashText(s: string): string {
  let h = 5381;
  for (let i = 0; i < s.length; i++) h = ((h << 5) + h + s.charCodeAt(i)) | 0;
  return (h >>> 0).toString(36) + ':' + s.length;
}

/** Auto-detect the drug ligand: the non-solvent HETATM residue with the most
 *  atoms. Returns its comp-id + chain, or empty strings if none. */
function detectLigand(pdbText: string): {compId: string; chain: string} {
  const counts = new Map<string, {compId: string; chain: string; n: number}>();
  for (const line of pdbText.split(/\r?\n/)) {
    if (!line.startsWith('HETATM')) continue;
    const resName = line.slice(17, 20).trim().toUpperCase();
    if (!resName || NON_LIGAND_HET.has(resName)) continue;
    const chain = line.slice(21, 22).trim();
    const resSeq = line.slice(22, 26).trim();
    const key = `${resName}|${chain}|${resSeq}`;
    const e = counts.get(key);
    if (e) e.n++;
    else counts.set(key, {compId: resName, chain, n: 1});
  }
  let best: {compId: string; chain: string; n: number} | null = null;
  for (const e of counts.values())
    if (!best || e.n > best.n) best = e;
  return best ? {compId: best.compId, chain: best.chain} : {compId: '', chain: ''};
}

/** Parse `REMARK   2 RESOLUTION.   2.00 ANGSTROMS.` → 2.0, else null. */
function parseResolution(pdbText: string): number | null {
  for (const line of pdbText.split(/\r?\n/)) {
    if (!line.startsWith('REMARK   2')) continue;
    const m = line.match(/RESOLUTION\.?\s+([0-9]+(?:\.[0-9]+)?)\s*ANGSTROM/i);
    if (m) return parseFloat(m[1]);
  }
  return null;
}

/** True if the text actually contains atomic coordinates. */
export function looksLikePdb(text: string): boolean {
  return /^(ATOM|HETATM)/m.test(text);
}

/** Register a dropped custom PDB file and return its registry entry. The
 *  returned `id` is what gets added to the textarea. */
export function registerLocalPdb(fileName: string, pdbText: string): LocalPdb {
  const lig = detectLigand(pdbText);
  const entry: LocalPdb = {
    id: makeLocalId(fileName),
    fileName,
    pdbText,
    resolution: parseResolution(pdbText),
    ligandCompId: lig.compId,
    ligandChain: lig.chain,
    contentHash: hashText(pdbText),
  };
  store.set(entry.id, entry);
  return entry;
}

/** Fingerprint contribution for the local ids present in `ids` (id + content
 *  hash). Lets the preview cache invalidate when a dropped file changes. */
export function localFingerprintFor(ids: string[]): string {
  return localIdsIn(ids)
    .map((id) => `${id.toUpperCase()}=${store.get(id.toUpperCase())!.contentHash}`)
    .sort().join(',');
}

/** Synthesized pdb_qc rows for the given local ids, matching the row shape
 *  `enrichPdbList` builds for RCSB entries. The drug ligand is the auto-detected
 *  HETATM residue; SMILES/MW are left empty (no reliable bond perception from
 *  coordinates), so the MW filter is effectively skipped for local structures. */
export function localQcRows(ids: string[]): Record<string, any>[] {
  return localIdsIn(ids).map((id) => {
    const e = store.get(id.toUpperCase())!;
    return {
      pdb_id: e.id,
      resolution: e.resolution,
      experimental_method: 'Local file',
      r_free: null,
      ligand_comp_id: e.ligandCompId,
      ligand_chain: e.ligandChain,
      ligand_formula_weight: null,
      ligand_smiles: '',
      ligand_rscc: null,
    };
  });
}
