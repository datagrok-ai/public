// Chem's RDKit substructure search vs Chem's crux engine, in Node, over a SMILES dataset and query sets.
//
//   node parity.mjs [molecules=20000] [sets=smiles,molblock,molblockH,self]
//
// sets:       smiles, molblock (sketcher-like molblocks of the smiles set), molblockH (with explicit hydrogens),
//             self (dataset molecules as queries), bdb-100 / bdb-1k / bdb-3488 (SMARTS sets from crux-bench)
// DATASET     SMILES csv (default: <reddata>/data/demo/chem/chembl/chembl-100k.csv)
// CRUX_BENCH  crux-bench checkout, for the bdb-* sets
// CRUX_WASM   a wasm-pack out dir to test instead of the build vendored in Chem
// PURE=1      leave out the RDKit verdict Chem gives the molecules crux cannot parse (shows crux's own results)
// VERBOSE=1   print every query, not only fallbacks and differences
import {readFileSync} from 'fs';
import {join} from 'path';
import {performance} from 'perf_hooks';
import {DATA_DIR, loadChemQueryCode, loadCrux, loadRdkit, readSmilesCsv} from './node-env.mjs';

const N = parseInt(process.argv[2] ?? '20000');
const SETS = (process.argv[3] ?? 'smiles,molblock,molblockH,self').split(',');
const DATASET = process.env.DATASET ?? join(DATA_DIR, 'demo', 'chem', 'chembl', 'chembl-100k.csv');
const PURE = process.env.PURE === '1';

const SMILES_QUERIES = ['c1ccccc1', 'c1ccncc1', 'C1CC1', 'C1CCCCC1', 'C1CCNCC1', 'C1CNCCN1', 'C1COCCN1',
  'c1ccc2[nH]ccc2c1', 'c1ccc2ncccc2c1', 'c1cn[nH]c1', 'c1c[nH]cn1', 'c1cscn1', 'c1ccoc1', 'c1ccsc1', 'C(=O)N', 'C(=O)O',
  'OC=O', 'S(=O)(=O)N', 'C#N', '[N+](=O)[O-]', 'C(F)(F)F', 'Cl', 'Br', 'F', 'I', 'N', 'O', 'S', 'P', 'CC', 'CCC', 'CCCC',
  'C=C', 'C#C', 'CO', 'CN', 'CS', 'NC(N)=O', 'c1ccc(cc1)-c1ccccc1', 'O=C1CCCN1', 'c1ccc2ccccc2c1', 'c1cc[nH]c1',
  'C[C@H](N)C(=O)O', 'N[C@@H](Cc1ccccc1)C(=O)O', '[NH3+]', '[O-]', 'C(=O)[O-]', '[Na+]', 'OCC(O)CO', 'c1ccc(O)cc1',
  'Nc1ccccc1', 'CC(C)C', 'C1=CC=CC=C1', 'C1=CN=CC=C1', 'O=C(Nc1ccccc1)c1ccccc1', 'c1ccc2c(c1)OCO2', 'C1CCC2CCCCC2C1',
  'c1ccc2c(c1)[nH]c1ccccc12', 'O=c1cc[nH]cc1', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'CC(=O)Nc1ccc(O)cc1', 'C/C=C/C',
  'F/C=C/F', 'C[C@@H](O)CC', 'c1ccc(-n2cccc2)cc1', 'C1CC2CCC1C2', 'C12CC3CC(CC(C3)C1)C2', '[Si]', '[B]', 'B', 'Se', '[se]',
  'OB(O)c1ccccc1', 'CS(C)=O', 'c1ccc2c(c1)CC(=O)N2', 'O=S(=O)(O)c1ccccc1', 'NC(=N)N', 'C=CC(=O)N', 'OO', 'NN', 'N=N',
  'N#N'];

const rdkit = await loadRdkit();
const {crux, dir: cruxDir} = await loadCrux();
const {getMolSafe, getQueryMolSafe, getCruxSmarts} = await loadChemQueryCode();
const mols = readSmilesCsv(DATASET, N);
console.log(`${mols.length} molecules from ${DATASET}\ncrux build: ${cruxDir}${PURE ? ' (PURE)' : ''}`);

// the RDKit search: canonical SMILES + pattern fingerprint per molecule (getFingerprints), then the query is
// matched against the canonical SMILES parsed without sanitization (searchWithPatternFps)
let t = performance.now();
const patternFps = new Array(mols.length).fill(null);
const targets = new Array(mols.length).fill(null);
const unsanitized = JSON.stringify({sanitize: false, removeHs: false, assignStereo: false});
for (let i = 0; i < mols.length; i++) {
  const mol = mols[i] ? getMolSafe(mols[i], {}, rdkit).mol : null;
  if (!mol)
    continue;
  try {
    patternFps[i] = mol.get_pattern_fp_as_uint8array();
    targets[i] = rdkit.get_mol(mol.get_smiles(), unsanitized);
  } catch (e) {
  } finally {
    mol.delete();
  }
}
console.log(`RDKit: prepared in ${(performance.now() - t).toFixed(0)} ms`);

function rdkitSearch(query) {
  const qmol = getQueryMolSafe(query, '', rdkit);
  if (!qmol)
    return null;
  try {
    const qfp = qmol.get_pattern_fp_as_uint8array();
    const hits = [];
    outer:
    for (let i = 0; i < targets.length; i++) {
      const fp = patternFps[i];
      if (!fp || !targets[i])
        continue;
      for (let j = 0; j < qfp.length; j++) {
        if ((fp[j] & qfp[j]) !== qfp[j])
          continue outer;
      }
      if (targets[i].get_substruct_match(qmol) !== '{}')
        hits.push(i);
    }
    return hits;
  } finally {
    qmol.delete();
  }
}

// Chem's crux engine: the index over the column (CXSMILES extensions cut), getCruxQuery, and the RDKit verdict
// for the molecules crux cannot parse
t = performance.now();
const builder = new crux.CollectionBuilder(true);
const unparsed = [];
for (let i = 0; i < mols.length; i++) {
  const extension = mols[i]?.indexOf(' |') ?? -1;
  if (builder.addMany([extension > 0 ? mols[i].substring(0, extension) : mols[i] ?? '']) === 0)
    unparsed.push(i);
}
const collection = builder.finish();
const unparsedMolecules = unparsed.filter((i) => {
  const mol = mols[i] ? getMolSafe(mols[i], {}, rdkit).mol : null;
  mol?.delete();
  return !!mol;
});
console.log(`crux: indexed in ${(performance.now() - t).toFixed(0)} ms, ${unparsed.length} rows not parsed, ` +
  `${unparsedMolecules.length} of them molecules for RDKit: ${unparsedMolecules.slice(0, 5).map((i) => mols[i]).join(' ')}`);
const validator = new crux.CollectionBuilder(true).finish();

function cruxQuery(query) {
  const qmol = getQueryMolSafe(query, '', rdkit);
  if (!qmol)
    return null;
  try {
    const smarts = getCruxSmarts(qmol);
    if (!smarts)
      return null;
    validator.substructureSearch(smarts, 0);
    return smarts;
  } catch (e) {
    return null;
  } finally {
    qmol.delete();
  }
}

function cruxSearch(query, smarts) {
  const hits = new Set(collection.substructureSearch(smarts, 0));
  if (!PURE && unparsedMolecules.length > 0) {
    const qmol = getQueryMolSafe(query, '', rdkit);
    for (const i of unparsedMolecules) {
      const mol = getMolSafe(mols[i], {}, rdkit).mol;
      if (mol.get_substruct_match(qmol) !== '{}')
        hits.add(i);
      mol.delete();
    }
    qmol.delete();
  }
  return Array.from(hits).sort((a, b) => a - b);
}

function molblock(smiles, explicitHydrogens) {
  const mol = rdkit.get_mol(smiles);
  if (!mol)
    return null;
  try {
    return explicitHydrogens ? mol.add_hs() : mol.get_molblock();
  } finally {
    mol.delete();
  }
}

const queries = [];
if (SETS.includes('smiles'))
  queries.push(...SMILES_QUERIES.map((q) => ({set: 'smiles', q, label: q})));
if (SETS.includes('molblock'))
  queries.push(...SMILES_QUERIES.map((s) => ({set: 'molblock', q: molblock(s, false), label: s})).filter((q) => q.q));
if (SETS.includes('molblockH'))
  queries.push(...SMILES_QUERIES.map((s) => ({set: 'molblockH', q: molblock(s, true), label: `${s} +H`})).filter((q) => q.q));
if (SETS.includes('self')) {
  for (let k = 0; k < 40; k++) {
    const i = Math.floor((k + 0.5) * mols.length / 40);
    if (targets[i])
      queries.push({set: 'self', q: mols[i], label: `row ${i}`});
  }
}
for (const set of SETS.filter((s) => s.startsWith('bdb-'))) {
  if (!process.env.CRUX_BENCH)
    throw new Error(`set ${set} needs CRUX_BENCH (a crux-bench checkout)`);
  const text = readFileSync(join(process.env.CRUX_BENCH, 'queries', set, 'queries.txt'), 'utf8');
  queries.push(...text.split(/\r?\n/).filter((l) => l.trim()).map((q) => ({set, q, label: q})));
}

let differences = 0;
let fallbacks = 0;
const total = {rdkit: 0, crux: 0};
for (const {set, q, label} of queries) {
  t = performance.now();
  const rdkitHits = rdkitSearch(q);
  const rdkitMs = performance.now() - t;
  total.rdkit += rdkitMs;
  const smarts = cruxQuery(q);
  t = performance.now();
  const cruxHits = smarts === null ? rdkitHits : cruxSearch(q, smarts);
  const cruxMs = performance.now() - t;
  total.crux += cruxMs;
  const inCrux = new Set(cruxHits ?? []);
  const inRdkit = new Set(rdkitHits ?? []);
  const onlyRdkit = (rdkitHits ?? []).filter((i) => !inCrux.has(i));
  const onlyCrux = (cruxHits ?? []).filter((i) => !inRdkit.has(i));
  const differs = onlyRdkit.length + onlyCrux.length > 0;
  differences += differs ? 1 : 0;
  fallbacks += smarts === null ? 1 : 0;
  if (!differs && smarts !== null && process.env.VERBOSE !== '1')
    continue;
  const flag = differs ? '!!' : smarts === null ? 'RD' : '  ';
  console.log(`${flag} [${set}] ${label.slice(0, 36).padEnd(36)} rdkit ${String(rdkitHits?.length).padStart(7)} ` +
    `crux ${String(cruxHits?.length).padStart(7)}  ${rdkitMs.toFixed(0).padStart(5)} / ${cruxMs.toFixed(1).padStart(6)} ms  ` +
    `${smarts ?? '(RDKit)'}`.slice(0, 250));
  if (differs) {
    console.log(`     only RDKit: ${onlyRdkit.slice(0, 3).map((i) => mols[i]).join(' ')}`);
    console.log(`     only crux:  ${onlyCrux.slice(0, 3).map((i) => mols[i]).join(' ')}`);
  }
}
console.log(`\n${queries.length} queries: ${differences} differ, ${fallbacks} run on RDKit (RD); ` +
  `search time RDKit ${total.rdkit.toFixed(0)} ms, crux ${total.crux.toFixed(0)} ms (one thread each)`);
process.exitCode = differences > 0 ? 1 : 0;
