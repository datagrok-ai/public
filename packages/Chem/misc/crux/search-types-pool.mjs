// The search types beside Contains over a pool of worker threads, as Chem's two engines split a column over their web
// workers, timed: each type's first search (Crux: the segments' keys; RDKit: the pool's pattern fingerprints and
// canonical SMILES) and a repeated one, with Chem's own code (compare100k.mjs's TypeScript loader).
//
//   node search-types-pool.mjs [rows=100000] [workers=14]
//
// RDKit mode: each worker holds a slice of the column (≤ 10,000 rows: its cache), RdKitServiceWorker.getFingerprints
// (Pattern, canonical SMILES) on the first search, then RdKitService's pattern screen and the worker's
// searchSubstructure on the rows it passes. Crux mode: crux/jobs/search-types.ts's job on segments laid out as
// CruxService lays out an index (1,000 rows per worker in the first round, doubling, the last round split evenly).
// Times: the engines' work over the pool, end to end in the threads (the browser's message passing and the main
// thread's merge are not in them). DATASET as compare100k.mjs.
import {Worker, isMainThread, parentPort, workerData} from 'worker_threads';
import {fileURLToPath} from 'url';
import {join, dirname} from 'path';
import {performance} from 'perf_hooks';

const HERE = dirname(fileURLToPath(import.meta.url));

if (isMainThread) {
  const rows = parseInt(process.argv[2] ?? '100000');
  const workers = parseInt(process.argv[3] ?? '14');
  const {molecules} = await import(join(HERE, 'compare100k.mjs'));
  const mols = molecules(rows);
  const types = ['Exact', 'Stereo agnostic', 'Included in', 'Not included in'];
  const queries = [mols[1234], mols[56789 % mols.length]];
  const pool = Array.from({length: workers}, () => new Worker(fileURLToPath(import.meta.url)));
  const call = (w, msg) => new Promise((resolve, reject) => {
    pool[w].once('message', (m) => m.error ? reject(new Error(m.error)) : resolve(m));
    pool[w].postMessage(msg);
  });
  // the segments of each worker
  const crux = Array.from({length: workers}, () => []);
  for (let start = 0, size = 1000; start < rows; size = Math.min(50000, size * 2)) {
    const round = Math.max(500, Math.min(size, Math.ceil((rows - start) / workers)));
    for (let w = 0; w < workers && start < rows; w++) {
      const end = Math.min(rows, start + round);
      crux[w].push([start, end]);
      start = end;
    }
  }
  const rdkit = Array.from({length: workers}, (_, w) => {
    const per = Math.ceil(rows / workers);
    return [[w * per, Math.min(rows, (w + 1) * per)]];
  });
  const lines = [];
  for (const engine of ['rdkit', 'crux']) {
    await Promise.all(pool.map((_, w) => call(w, {op: 'init', engine})));
    for (const type of types) {
      // a fresh column per type: each starts without keys / fingerprints
      await Promise.all(pool.map((_, w) => call(w, {op: 'reset'})));
      const times = [];
      for (const query of queries) {
        const t0 = performance.now();
        const parts = await Promise.all(pool.map((_, w) => call(w, {op: 'search', type, query,
          segments: (engine === 'crux' ? crux : rdkit)[w].map(([s, e]) => ({start: s, molecules: mols.slice(s, e)}))})));
        const ms = performance.now() - t0;
        times.push(`${Math.round(ms)} ms (${parts.reduce((t, p) => t + p.hits, 0)} hits)`);
      }
      lines.push(`${engine} ${type}: first ${times[0]}, repeated ${times[1]}`);
      console.log(lines[lines.length - 1]);
    }
  }
  await Promise.all(pool.map((w) => w.terminate()));
  console.log(`${rows} rows, ${workers} workers`);
} else {
  const h = await import(join(HERE, 'compare100k.mjs'));
  let state = null;
  parentPort.on('message', async (msg) => {
    try {
      if (msg.op === 'init') {
        const {module} = await h.engineModules(msg.engine);
        const jobs = await h.chem('crux/jobs/search-types.ts');
        const {RdKitServiceWorker} = await h.chem('rdkit-service/rdkit-service-worker.ts');
        const {getQueryMolSafe} = await h.chem('utils/mol-creation_rdkit.ts');
        state = {engine: msg.engine, module, jobs, RdKitServiceWorker, getQueryMolSafe, segments: new Map()};
        parentPort.postMessage({ok: true});
      } else if (msg.op === 'reset') {
        state.segments.clear();
        parentPort.postMessage({ok: true});
      } else if (msg.op === 'search') {
        let hits = 0;
        for (const {start, molecules} of msg.segments) {
          if (state.engine === 'crux') {
            // the job of crux.worker.ts: the segment's molecules on its first search
            const known = state.segments.has(start);
            state.segments.set(start, true);
            const res = state.jobs.SEARCH_TYPE_JOBS['searchTypes::search']({module: state.module}, String(start),
              known ? null : molecules, msg.query, '', msg.type);
            hits += res.hits.length;
          } else {
            let seg = state.segments.get(start);
            if (!seg) {
              const worker = new state.RdKitServiceWorker(state.module, '');
              seg = {worker, fps: await worker.getFingerprints('Pattern', molecules, true)};
              state.segments.set(start, seg);
            }
            const qmol = state.getQueryMolSafe(msg.query, '', state.module);
            const qfp = qmol.get_pattern_fp_as_uint8array();
            qmol.delete();
            const not = msg.type === 'Not included in';
            const passed = [];
            for (let i = 0; i < seg.fps.fps.length; i++) {
              const fp = seg.fps.fps[i];
              if (not) {
                passed.push(i);
                continue;
              }
              if (!fp)
                continue;
              let ok = true;
              for (let j = 0; j < 256 && ok; j++)
                ok = (fp[j] & qfp[j]) === (msg.type === 'Included in' ? fp[j] : qfp[j]);
              if (ok)
                passed.push(i);
            }
            const bits = new Uint32Array(await seg.worker.searchSubstructure(msg.query, '',
              passed.map((i) => seg.fps.smiles[i]), msg.type));
            passed.forEach((_, k) => hits += bits[k >>> 5] & (1 << (k & 31)) ? 1 : 0);
          }
        }
        parentPort.postMessage({hits});
      }
    } catch (e) {
      parentPort.postMessage({error: String(e?.stack ?? e)});
    }
  });
}
