// crux memory per indexed molecule, index build speed and the molecules crux cannot parse, over a SMILES csv.
//
//   node mem.mjs <csv> [molecules] [segment size=33000]
//
// CRUX_WASM  a wasm-pack out dir to test instead of the build vendored in Chem
import {performance} from 'perf_hooks';
import {loadCrux, readSmilesCsv} from './node-env.mjs';

const [csv, limit, segmentArg] = process.argv.slice(2);
if (!csv)
  throw new Error('usage: node mem.mjs <csv> [molecules] [segment size]');
const SEGMENT = parseInt(segmentArg ?? '33000');
const {crux, wasm, dir} = await loadCrux();
const mols = readSmilesCsv(csv, limit ? parseInt(limit) : Infinity);

const before = wasm.memory.buffer.byteLength;
const collections = [];
const failed = [];
let t = performance.now();
// segments like Chem's largest ones; wasm memory never shrinks, so this also counts the build high-water mark
for (let start = 0; start < mols.length; start += SEGMENT) {
  const builder = new crux.CollectionBuilder(true);
  builder.setLenient?.(true);
  mols.slice(start, start + SEGMENT).forEach((smiles, i) => {
    if (builder.addMany([smiles ?? '']) === 0)
      failed.push(start + i);
  });
  collections.push(builder.finish());
}
const buildMs = performance.now() - t;
const bytes = wasm.memory.buffer.byteLength - before;
console.log(`crux build: ${dir}`);
console.log(`${mols.length} molecules in ${collections.length} segments: index built in ${buildMs.toFixed(0)} ms ` +
  `(${(1000 * buildMs / mols.length).toFixed(1)} us/mol), wasm memory +${(bytes / 1e6).toFixed(1)} MB ` +
  `(${(bytes / mols.length).toFixed(0)} B/mol)`);
console.log(`${failed.length} molecules not parsed: ${failed.slice(0, 10).map((i) => mols[i]).join(' ')}`);

for (const smarts of ['[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1', '[#6]1-[#6]-[#6]-1', '[#7]', '[#6](=[#8])-[#7]']) {
  t = performance.now();
  let hits = 0;
  for (const collection of collections)
    hits += collection.substructureSearch(smarts, 0).length;
  console.log(`${smarts.padEnd(36)} ${String(hits).padStart(8)} hits in ${(performance.now() - t).toFixed(0)} ms`);
}
