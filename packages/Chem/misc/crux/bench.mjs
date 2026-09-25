// Chem substructure search in the browser, on a running stand with Chem published: RDKit vs Crux engine, through
// grok.chem.searchSubstructure (timings and hit sets) and through the substructure filter (time to the first
// result, time to the end, filter updates).
//
//   node bench.mjs <csv> [modes=api,filter] [engines=RDKit,Crux]
//
// STAND_API  the server's API URL, for the dev-key login (default http://localhost:8082)
// STAND_WEB  the client on the same origin as the API: package web workers cannot start cross-origin, so not the
//            pub serve port (default http://localhost:8888)
// DEV_KEY    developer key (default admin)
// It sets the user's Substructure Search Engine package property and leaves it on RDKit.
import {writeFileSync} from 'fs';
import {loadRdkit, requireFromChem} from './node-env.mjs';

const [csvPath, modesArg, enginesArg] = process.argv.slice(2);
if (!csvPath)
  throw new Error('usage: node bench.mjs <csv> [api,filter] [RDKit,Crux]');
const modes = (modesArg ?? 'api,filter').split(',');
const engines = (enginesArg ?? 'RDKit,Crux').split(',');
const API = process.env.STAND_API ?? 'http://localhost:8082';
const WEB = process.env.STAND_WEB ?? 'http://localhost:8888';

const rdkit = await loadRdkit();
const molblock = (smiles) => {
  const mol = rdkit.get_mol(smiles);
  try {
    return mol.get_molblock();
  } finally {
    mol.delete();
  }
};
const QUERIES = [
  {label: 'benzene (smiles)', q: 'c1ccccc1'},
  {label: 'pyridine (smiles)', q: 'c1ccncc1'},
  {label: 'cyclopropane (smiles)', q: 'C1CC1'},
  {label: 'amide (smiles)', q: 'C(=O)N'},
  {label: 'indole [nH] (smiles)', q: 'c1ccc2[nH]ccc2c1'},
  {label: 'piperazine (molblock)', q: molblock('C1CNCCN1')},
  {label: 'quinoline (molblock)', q: molblock('c1ccc2ncccc2c1')},
  {label: 'sulfonamide (molblock)', q: molblock('NS(=O)(=O)c1ccccc1')},
  {label: 'alanine stereo (molblock)', q: molblock('C[C@H](N)C(=O)O')},
  {label: 'CF3 (molblock)', q: molblock('C(F)(F)F')},
  {label: 'carboxylic acid (smarts)', q: '[CX3](=O)[OX2H1]'},
  {label: 'amine not amide (smarts)', q: '[NX3;H2,H1;!$(NC=O)]'},
  {label: 'halogen (smarts)', q: '[F,Cl,Br,I]'},
  {label: '[OH] radical (RDKit either way)', q: '[OH]'},
];
const FILTER_QUERIES = [
  {label: 'benzene', q: molblock('c1ccccc1')},
  {label: 'piperazine', q: molblock('C1CNCCN1')},
  {label: 'thiazole', q: molblock('c1cscn1')},
];

const login = await fetch(`${API}/users/login/dev`, {method: 'POST',
  headers: {Authorization: `Dev ${process.env.DEV_KEY ?? 'admin'}`}});
const {token} = await login.json();
const {chromium} = requireFromChem('@playwright/test');
const browser = await chromium.launch({headless: true});
const page = await browser.newPage({viewport: {width: 1600, height: 1000}});
page.setDefaultTimeout(1800000);
page.on('console', (m) => {
  if (m.type() === 'error' || /crux/i.test(m.text()))
    console.log(`  [browser ${m.type()}] ${m.text().slice(0, 300)}`);
});
await page.route('**/__crux-bench/data.csv', (route) => route.fulfill({path: csvPath, contentType: 'text/csv'}));
await page.goto(`${WEB}/login.html?token=${encodeURIComponent(token)}`);
await page.waitForFunction(() => {
  try {
    return !!(window.grok?.shell?.user && DG.Func.find({package: 'Chem', name: 'searchSubstructure'}).length);
  } catch (e) {
    return false;
  }
}, null, {timeout: 300000, polling: 1000});

await page.evaluate(async () => {
  window.__csv = await (await fetch('/__crux-bench/data.csv')).text();
  window.__setEngine = async (engine) => {
    const pkg = DG.Func.find({package: 'Chem', name: 'searchSubstructure'})[0].package;
    await pkg.setSettings({SubstructureSearchEngine: engine}, grok.shell.user.group);
  };
  // the molecules: the first column with smiles in its name, else the detected one
  window.__newDf = async (name) => {
    const df = DG.DataFrame.fromCsv(window.__csv);
    df.name = name;
    await grok.data.detectSemanticTypes(df);
    const col = df.columns.toList().find((c) => /smiles/i.test(c.name)) ?? df.columns.bySemType(DG.SEMTYPE.MOLECULE);
    col.semType = DG.SEMTYPE.MOLECULE;
    return {df, col};
  };
});

const results = {api: {}, filter: {}};
try {
  if (modes.includes('api')) {
    for (const engine of engines) {
      results.api[engine] = JSON.parse(await page.evaluate(async ({engine, queries}) => {
        const hash = (buf) => {
          let h = 0x811c9dc5;
          for (const x of buf)
            h = Math.imul(h ^ x, 0x01000193);
          return (h >>> 0).toString(16);
        };
        await window.__setEngine(engine);
        const {df, col} = await window.__newDf(`bench_api_${engine}`);
        let t = performance.now();
        await grok.chem.searchSubstructure(col, 'c1ccccc1');
        const out = [{label: 'first search (fingerprints / index)', ms: performance.now() - t}];
        for (const {label, q} of queries) {
          t = performance.now();
          const bs = await grok.chem.searchSubstructure(col, q);
          out.push({label, ms: performance.now() - t, count: bs.trueCount, hash: hash(bs.getBuffer())});
        }
        return JSON.stringify({rows: df.rowCount, out});
      }, {engine, queries: QUERIES}));
    }
    for (const engine of engines) {
      console.log(`\n[api] ${engine}, ${results.api[engine].rows} rows`);
      for (const r of results.api[engine].out)
        console.log(`  ${r.label.padEnd(40)} ${r.ms.toFixed(0).padStart(7)} ms  ${r.count ?? ''}`);
    }
    if (engines.length === 2) {
      const [a, b] = engines.map((e) => results.api[e].out);
      console.log(`\n[api] ${engines.join(' vs ')}`);
      for (let i = 1; i < a.length; i++) {
        const same = a[i].hash === b[i].hash && a[i].count === b[i].count;
        console.log(`  ${same ? 'same' : 'DIFF'}  ${a[i].label.padEnd(40)} ${String(a[i].count).padStart(8)} vs ` +
          `${String(b[i].count).padStart(8)}   ${a[i].ms.toFixed(0).padStart(6)} ms vs ${b[i].ms.toFixed(0).padStart(6)} ms`);
      }
    }
  }

  if (modes.includes('filter')) {
    for (const engine of engines) {
      results.filter[engine] = JSON.parse(await page.evaluate(async ({engine, queries}) => {
        await window.__setEngine(engine);
        // the filter is created synchronously, the package must be loaded before
        await grok.functions.call('Chem:substructureFilter');
        const {df, col} = await window.__newDf(`bench_filter_${engine}`);
        const tv = grok.shell.addTableView(df);
        const colName = col.name;
        const runs = [];
        for (const {label, q} of queries) {
          // the search starts when the filter announces its new query (it applies a state after a second)
          const run = {label, start: 0, progress: [], filterChanges: [], done: 0};
          let armed = false;
          const subs = [
            grok.events.onCustomEvent('chem-substructure-filter').subscribe((s) => {
              if (armed && !run.start && s.molblock)
                run.start = performance.now();
            }),
            grok.events.onCustomEvent(`substructure_search_progress-${df.name}-${colName}`).subscribe((p) => {
              if (run.start)
                run.progress.push({t: performance.now() - run.start, p: Number(p)});
            }),
            df.onFilterChanged.subscribe(() => {
              if (run.start)
                run.filterChanges.push({t: performance.now() - run.start, count: df.filter.trueCount});
            }),
            grok.events.onCustomEvent(`terminate_substructure_search-${df.name}-${colName}`).subscribe((key) => {
              if (run.start && key && !run.done)
                run.done = performance.now() - run.start;
            }),
          ];
          armed = true;
          tv.getFiltersGroup({createDefaultFilters: false}).updateOrAdd({
            type: DG.FILTER_TYPE.SUBSTRUCTURE, column: colName, columnName: colName, molBlock: q}, false);
          const t0 = performance.now();
          while ((!run.done || performance.now() - run.start - run.done < 1500) && performance.now() - t0 < 900000)
            await new Promise((r) => setTimeout(r, 50));
          subs.forEach((s) => s.unsubscribe());
          run.finalCount = df.filter.trueCount;
          runs.push(run);
        }
        grok.shell.closeAll();
        return JSON.stringify({rows: df.rowCount, runs});
      }, {engine, queries: FILTER_QUERIES}));
      const {rows, runs} = results.filter[engine];
      console.log(`\n[filter] ${engine}, ${rows} rows`);
      for (const r of runs) {
        const first = r.filterChanges.find((c) => c.count > 0 && c.count < rows) ?? r.filterChanges[0];
        console.log(`  ${r.label.padEnd(12)} first result ${first ? first.t.toFixed(0) : '-'} ms, ` +
          `done ${r.done.toFixed(0)} ms, ${r.filterChanges.length} filter updates, final ${r.finalCount}`);
      }
    }
  }
} finally {
  await page.evaluate(() => window.__setEngine('RDKit'));
  await browser.close();
}
const out = process.env.OUT ?? 'crux-bench-results.json';
writeFileSync(out, JSON.stringify(results, null, 1));
console.log(`\nresults: ${out}`);
