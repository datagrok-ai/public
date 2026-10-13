// Where crux and Chem's RDKit disagree on the molecules of a dataset, molecule by molecule:
//   A  rows RDKit reads the way Chem parses a column (getMolSafe) and crux cannot parse,
//   B  rows crux parses and RDKit cannot,
//   C  rows both parse but crux reads differently: crux does not find the molecule with RDKit's own canonical SMILES
//      as the query (translated by getCruxSmarts, as Chem translates queries) while RDKit does.
//
//   node census.mjs <csv|smi> [molecules]
//
// THREADS    worker threads (default: cores - 2)
// OUT        a TSV of every disagreement
// CRUX_WASM  a wasm-pack out dir to test instead of the build vendored in Chem
import {writeFileSync} from 'fs';
import {availableParallelism} from 'os';
import {performance} from 'perf_hooks';
import {Worker, isMainThread, parentPort, workerData} from 'worker_threads';
import {loadChemQueryCode, loadCrux, loadRdkit, readSmilesCsv} from './node-env.mjs';

// what a disagreement is about, from the input text (first match wins)
const CLASSES = [
  ['element beyond Rn', /\[\d*(Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og)[^a-z]/],
  ['aromatic N-oxide n(=O)', /n\(=O\)|O=n|n=O/],
  ['nitro/N-oxide N(=O)=', /N\(=O\)=|N\(=O\)\(=O\)|O=N\(=O\)|O=N=|=N\(=O\)/],
  ['aromatic p/as', /(^|[^A-Z[])(p|as)[\d(]|\[(p|as)[^a-z]/],
  ['wildcard *', /\*/],
  ['explicit [H] atom', /\[\d*H[+-]?\]/],
  ['quadruple bond $', /\$/],
  ['aromatic bond :', /[a-zA-Z\]]:[a-zA-Z[]/],
  ['atomic number [#n]', /\[\d*#\d/],
];

function classify(smiles) {
  return CLASSES.find(([, re]) => re.test(smiles))?.[0] ?? 'other';
}

async function censusRows(rows) {
  const rdkit = await loadRdkit();
  const {crux} = await loadCrux();
  const {getMolSafe, getQueryMolSafe, getCruxSmarts} = await loadChemQueryCode();
  const unsanitized = JSON.stringify({sanitize: false, removeHs: false, assignStereo: false});
  const counts = {};
  const count = (key) => counts[key] = (counts[key] ?? 0) + 1;
  const records = [];
  for (const [row, smiles] of rows) {
    if (!smiles) {
      count('empty');
      continue;
    }
    const {mol, kekulize, isQMol} = getMolSafe(smiles, {}, rdkit);
    const stage = !mol ? 'none' : isQMol ? 'qmol' : kekulize ? 'default' : 'kekulize:false';
    count(`rdkit ${stage}`);
    // crux gets the column value the way Chem's crux index reads it: CXSMILES extension cut
    const extension = smiles.indexOf(' |');
    const builder = new crux.CollectionBuilder(true);
    builder.setLenient?.(true);
    const cruxOk = builder.addMany([extension > 0 ? smiles.substring(0, extension) : smiles]) === 1;
    const collection = builder.finish();
    count(cruxOk ? 'crux parsed' : 'crux not parsed');
    try {
      if (mol && !cruxOk) {
        records.push({row, smiles, kind: 'A', stage, cls: classify(smiles)});
        continue;
      }
      if (!mol) {
        if (cruxOk)
          records.push({row, smiles, kind: 'B', stage, cls: classify(smiles)});
        continue;
      }
      if (isQMol) {
        count('self: not checked (RDKit reads the row as SMARTS)');
        continue;
      }
      const canonical = mol.get_smiles();
      const qmol = getQueryMolSafe(canonical, '', rdkit);
      const smarts = qmol ? getCruxSmarts(qmol) : null;
      let cruxHit = null;
      try {
        cruxHit = smarts === null ? null : collection.substructureSearch(smarts, 0).length === 1;
      } catch {
        // crux cannot parse the query: Chem runs such a search on RDKit
      }
      if (cruxHit === null) {
        count('self: not checked (query runs on RDKit)');
        qmol?.delete();
        continue;
      }
      const target = rdkit.get_mol(canonical, unsanitized);
      const rdkitHit = !!target && target.get_substruct_match(qmol) !== '{}';
      target?.delete();
      qmol.delete();
      if (!rdkitHit)
        count('self: not checked (RDKit misses itself)');
      else if (!cruxHit)
        records.push({row, smiles, kind: 'C', stage, cls: classify(smiles), canonical, smarts});
      else
        count('self: crux finds it');
    } finally {
      mol?.delete();
      collection.free();
    }
  }
  return {counts, records};
}

if (!isMainThread)
  parentPort.postMessage(await censusRows(workerData));
else {
  const [path, limit] = process.argv.slice(2);
  if (!path)
    throw new Error('usage: node census.mjs <csv|smi> [molecules]');
  const mols = readSmilesCsv(path, limit ? parseInt(limit) : Infinity);
  const threads = Math.max(1, Math.min(mols.length, parseInt(process.env.THREADS ?? '0') || availableParallelism() - 2));
  const t = performance.now();
  const results = await Promise.all(Array.from({length: threads}, (_, k) => new Promise((resolve, reject) => {
    const rows = [];
    for (let i = k; i < mols.length; i += threads)
      rows.push([i, mols[i]]);
    const worker = new Worker(new URL(import.meta.url), {workerData: rows});
    worker.once('message', resolve);
    worker.once('error', reject);
  })));
  const counts = {};
  for (const r of results) {
    for (const [key, n] of Object.entries(r.counts))
      counts[key] = (counts[key] ?? 0) + n;
  }
  const records = results.flatMap((r) => r.records).sort((a, b) => a.row - b.row);
  const {dir} = await loadCrux();
  const n = (key) => counts[key] ?? 0;
  console.log(`${mols.length} rows from ${path} (${n('empty')} empty), crux build ${dir}, ` +
    `${threads} threads, ${((performance.now() - t) / 1000).toFixed(0)} s`);
  console.log(`RDKit (getMolSafe): default ${n('rdkit default')}, kekulize:false ${n('rdkit kekulize:false')}, ` +
    `as SMARTS ${n('rdkit qmol')}, not parsed ${n('rdkit none')}`);
  console.log(`crux: parsed ${n('crux parsed')}, not parsed ${n('crux not parsed')}`);
  const titles = {A: 'crux cannot parse, RDKit can', B: 'crux parses, RDKit cannot',
    C: 'crux reads differently (misses the molecule by its RDKit SMILES)'};
  for (const kind of ['A', 'B', 'C']) {
    const own = records.filter((r) => r.kind === kind);
    console.log(`${kind}  ${titles[kind]}: ${own.length}`);
    const byClass = new Map();
    for (const r of own)
      byClass.set(r.cls, [...(byClass.get(r.cls) ?? []), r]);
    for (const [cls, list] of [...byClass].sort((a, b) => b[1].length - a[1].length)) {
      const stages = [...new Set(list.map((r) => r.stage))].join('/');
      console.log(`     ${String(list.length).padStart(6)}  ${cls} (RDKit ${stages}): ` +
        list.slice(0, 3).map((r) => r.smiles).join('  ').slice(0, 300));
    }
  }
  console.log(`self query: crux finds ${n('self: crux finds it')}; not checked: ` +
    `${n('self: not checked (query runs on RDKit)')} run on RDKit, ` +
    `${n('self: not checked (RDKit reads the row as SMARTS)')} SMARTS rows, ` +
    `${n('self: not checked (RDKit misses itself)')} RDKit misses itself`);
  if (process.env.OUT) {
    writeFileSync(process.env.OUT, ['row\tkind\tclass\trdkit\tsmiles\trdkit_canonical\tcrux_query'].concat(
      records.map((r) => [r.row, r.kind, r.cls, r.stage, r.smiles, r.canonical ?? '', r.smarts ?? ''].join('\t')))
      .join('\n') + '\n');
  }
}
