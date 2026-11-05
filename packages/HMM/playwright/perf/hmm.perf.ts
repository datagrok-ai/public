/* Performance at scale: 100,000 antibodies generated from Bio's antibodies sample
 * (random point mutations and small indels of random parents, seeded), then the
 * package's features and Bio's WebLogo, timed in the browser. Runs only on request:
 *   HMM_PERF=1 grok test --host localhost --skip-build --skip-puppeteer --category perf
 * Results: playwright/test-output/perf/perf.json (and the console). HMM_PERF_ROWS
 * sets the row count (default 100000). */
import {test, Page} from '@playwright/test';
import * as fs from 'fs';
import {
  loginToDatagrok, specTestOptions, softStep, stepErrors,
} from '@datagrok-libraries/test/src/playwright/spec-login';

test.use(specTestOptions);

const ROWS = Number(process.env.HMM_PERF_ROWS ?? 100_000);
const OUT = 'playwright/test-output/perf';
const results: {step: string; ms: number; note?: string; heapMB?: number}[] = [];

async function record(page: Page, step: string, ms: number, note?: string): Promise<void> {
  const heapMB = await page.evaluate(() => Math.round(((performance as any).memory?.usedJSHeapSize ?? 0) / 2 ** 20));
  results.push({step, ms: Math.round(ms), note, heapMB});
  console.log(`${step.padEnd(52)} ${String(Math.round(ms)).padStart(8)} ms  heap ${heapMB} MB  ${note ?? ''}`);
}

/** Time `body` in the page; it returns an optional note. */
async function timed(page: Page, step: string, body: string, arg: any = null): Promise<void> {
  const {ms, note} = await page.evaluate(async ({body, arg}) => {
    const fn = new Function('arg', `return (async () => { ${body} })();`);
    const t0 = performance.now();
    const note = await fn(arg);
    return {ms: performance.now() - t0, note: note == null ? undefined : String(note)};
  }, {body, arg});
  await record(page, step, ms, note);
}

/** Wait (in page code) until a WebLogo has drawn: positions computed, no render pending. */
const WEBLOGO_SETTLED = `
  const settled = async (v, min = 1) => {
    for (let k = 0; k < 6000; k++) {
      await new Promise((r) => setTimeout(r, 10));
      if ((v.positions?.length ?? 0) >= min && !v.isRenderPending && v.getWidgetStatus().values['positions shown'] > 0)
        return;
    }
    throw new Error('WebLogo did not settle');
  };`;

test('HMM performance at scale', async ({page}) => {
  test.setTimeout(3_600_000);
  stepErrors.length = 0;
  fs.mkdirSync(OUT, {recursive: true});
  await loginToDatagrok(page);
  await page.evaluate(() => {
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  });

  await softStep(`generate ${ROWS} antibodies and open them`, async () => {
    await timed(page, `generate ${ROWS} rows (in page)`, `
      const src = await grok.dapi.files.readCsv('System:AppData/Bio/samples/antibodies.csv');
      const hc = src.col('AntibodyHC').toList(), lc = src.col('AntibodyLC').toList(), y = src.col('Y').toList();
      let seed = 20261011;
      const rnd = () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
        t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
      const AA = 'ACDEFGHIKLMNPQRSTVWY';
      const aa = () => AA[Math.floor(rnd() * 20)];
      const mutate = (s) => {
        let out = '';
        for (const c of s) out += rnd() < 0.03 ? aa() : c;
        for (let k = 0; k < 2; k++) {
          const at = Math.floor(rnd() * out.length);
          if (rnd() < 0.2) out = out.slice(0, at) + aa().repeat(1 + Math.floor(rnd() * 3)) + out.slice(at);
          else if (rnd() < 0.2) out = out.slice(0, at) + out.slice(at + 1 + Math.floor(rnd() * 3));
        }
        return out;
      };
      const n = arg, ids = new Array(n), hcs = new Array(n), lcs = new Array(n), ys = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        const p = Math.floor(rnd() * hc.length), q = Math.floor(rnd() * lc.length);
        ids[i] = 'syn' + i; hcs[i] = mutate(hc[p]); lcs[i] = mutate(lc[q]); ys[i] = (y[p] ?? 0) + rnd() - 0.5;
      }
      window.__perfDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('Antibody_ID', ids),
        DG.Column.fromStrings('AntibodyHC', hcs), DG.Column.fromStrings('AntibodyLC', lcs),
        DG.Column.fromFloat32Array('Y', ys)]);
      window.__perfDf.name = 'antibodies100k';
      return Math.round(hcs.reduce((s, x) => s + x.length, 0) / n) + ' residues per heavy chain';`, ROWS);
    await timed(page, 'detect semantic types', `
      await grok.data.detectSemanticTypes(window.__perfDf);
      return window.__perfDf.col('AntibodyHC').semType + ' ' + window.__perfDf.col('AntibodyHC').meta.units;`);
    await timed(page, 'open table view (grid drawn)', `
      grok.shell.addTableView(window.__perfDf);
      for (let k = 0; k < 3000 && !document.querySelector('[name="viewer-Grid"] canvas'); k++)
        await new Promise((r) => setTimeout(r, 10));
      await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));`);
  });

  await softStep('grid scrolling', async () => {
    await timed(page, 'grid: scroll to 10 places, each drawn', `
      const grid = grok.shell.tv.grid;
      const slow = [];
      for (let k = 1; k <= 10; k++) {
        const t0 = performance.now();
        grid.scrollToCell('AntibodyHC', Math.floor(window.__perfDf.rowCount * k / 11));
        await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
        slow.push(Math.round(performance.now() - t0));
      }
      grid.scrollToCell('AntibodyHC', 0);
      return 'per scroll ms: ' + slow.join(' ');`);
  });

  await softStep('WebLogo of the raw heavy chains', async () => {
    await timed(page, 'WebLogo: add on AntibodyHC (100K rows)', `${WEBLOGO_SETTLED}
      const v = await window.__perfDf.plot.fromType('WebLogo', {sequenceColumnName: 'AntibodyHC'});
      grok.shell.tv.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'WebLogo HC', 0.3);
      window.__wl = v;
      await settled(v);
      return v.positions.length + ' positions';`);
    await timed(page, 'WebLogo: filter to half the rows', `${WEBLOGO_SETTLED}
      window.__perfDf.filter.init((i) => i % 2 === 0);
      await new Promise((r) => setTimeout(r, 50));
      await settled(window.__wl);`);
    await timed(page, 'WebLogo: clear the filter', `${WEBLOGO_SETTLED}
      window.__perfDf.filter.setAll(true);
      await new Promise((r) => setTimeout(r, 50));
      await settled(window.__wl);`);
    await timed(page, 'WebLogo: value column Y, average', `${WEBLOGO_SETTLED}
      window.__wl.setOptions({valueColumnName: 'Y', valueAggrType: 'avg'});
      await new Promise((r) => setTimeout(r, 50));
      await settled(window.__wl);`);
    await page.evaluate(() => (window as any).__wl.close());
  });

  await softStep('align to V-set (hmmalign)', async () => {
    await timed(page, 'hmmalign: 100K heavy chains to Pfam V-set', `
      const df = window.__perfDf;
      const col = await grok.functions.call('HMM:alignToProfile',
        {sequences: df.col('AntibodyHC'), profile: 'Pfam family or accession', family: 'V-set', trim: true});
      col.name = 'HC (V-set)';
      df.columns.add(col);
      return col.get(0).length + ' columns';`);
    await timed(page, 'WebLogo: data logo of the aligned column', `${WEBLOGO_SETTLED}
      const v = await window.__perfDf.plot.fromType('WebLogo', {sequenceColumnName: 'HC (V-set)'});
      grok.shell.tv.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'WebLogo aligned', 0.3);
      await settled(v);
      v.close();`);
    await timed(page, 'WebLogo: V-set profile on the aligned column', `${WEBLOGO_SETTLED}
      const v = await window.__perfDf.plot.fromType('WebLogo', {profile: 'V-set', sequenceColumnName: 'HC (V-set)'});
      grok.shell.tv.dockManager.dock(v, DG.DOCK_TYPE.DOWN, null, 'V-set logo', 0.3);
      await settled(v, 109);
      window.__wlp = v;
      return v.positionWidth + ' px per position';`);
    await timed(page, 'WebLogo profile: click a letter (select rows)', `
      const v = window.__wlp;
      const status = v.getWidgetStatus();
      const [name, a] = Object.entries(status.hitAreas).filter(([k]) => k.startsWith('monomer '))
        .sort(([, p], [, q]) => q.height - p.height)[0];
      const r = status.parts.canvas.getBoundingClientRect();
      const at = {clientX: r.left + a.x + a.width / 2, clientY: r.top + a.y + a.height / 2, bubbles: true};
      status.parts.canvas.dispatchEvent(new MouseEvent('mousedown', at));
      return name + ': ' + window.__perfDf.selection.trueCount + ' rows selected';`);
    await page.evaluate(() => {
      (window as any).__wlp.close();
      (window as any).__perfDf.selection.setAll(false);
    });
  });

  await softStep('searches', async () => {
    await timed(page, 'Profile HMM Search: V-set over 100K (function)', `
      const df = window.__perfDf;
      const r = await grok.functions.call('HMM:searchWithHmm',
        {table: df, sequence: df.col('AntibodyHC'), model: 'PF07686', evalue: 10, gathering: false});
      return r.columns.byIndex(0).stats.valueCount + ' reported';`);
    await timed(page, 'Homology search (phmmer): row 1 over 100K', `
      const df = window.__perfDf;
      const r = await grok.functions.call('HMM:homologySearch',
        {table: df, sequence: df.col('AntibodyHC'), query: df.col('AntibodyHC').get(0), evalue: 10});
      return r.columns.byIndex(0).stats.valueCount + ' reported';`);
    await timed(page, 'Find Domains: built-in library, GA, over 100K (function)', `
      const df = window.__perfDf;
      const r = await grok.functions.call('HMM:findDomains',
        {table: df, sequence: df.col('AntibodyHC'), library: 'builtin', evalue: 0.01, gathering: true});
      return r.rowCount + ' domains';`);
  });

  await softStep('Find Domains through the dialog (annotations written)', async () => {
    const t0 = Date.now();
    await page.evaluate(() => grok.functions.call('HMM:findDomainsDialog'));
    const dlg = page.locator('[name="dialog-Find-Domains-(HMMER)"]');
    await dlg.waitFor({timeout: 60_000});
    const t1 = Date.now();
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC domains'), null,
      {timeout: 1_800_000, polling: 200});
    await record(page, 'Find Domains dialog: OK to annotated column', Date.now() - t1,
      `dialog opened in ${t1 - t0} ms`);
    await timed(page, 'grid: scroll with domain annotations', `
      const grid = grok.shell.tv.grid;
      const slow = [];
      for (let k = 1; k <= 5; k++) {
        const t0 = performance.now();
        grid.scrollToCell('AntibodyHC', Math.floor(window.__perfDf.rowCount * k / 6));
        await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
        slow.push(Math.round(performance.now() - t0));
      }
      return 'per scroll ms: ' + slow.join(' ');`);
  });

  await softStep('build a profile', async () => {
    await timed(page, 'hmmbuild: 1,000 aligned rows', `
      const df = window.__perfDf;
      df.selection.init((i) => i < 1000);
      const text = await grok.functions.call('HMM:buildProfile',
        {table: df, alignment: df.col('HC (V-set)'), name: 'p1k', selectedOnly: true});
      return /^LENG\\s+(\\d+)/m.exec(text)[1] + ' positions';`);
    await timed(page, 'hmmbuild: all 100K aligned rows', `
      const df = window.__perfDf;
      df.selection.setAll(false);
      const text = await grok.functions.call('HMM:buildProfile',
        {table: df, alignment: df.col('HC (V-set)'), name: 'pall', selectedOnly: false});
      return /^LENG\\s+(\\d+)/m.exec(text)[1] + ' positions';`);
  });

  await softStep('ANARCI numbering', async () => {
    await timed(page, 'ANARCI engine: 100K heavy chains (IMGT)', `
      const df = window.__perfDf;
      const r = await grok.functions.call('HMM:anarciNumbering',
        {df, seqCol: df.col('AntibodyHC'), scheme: 'imgt'});
      return r.col('chain').stats.valueCount + ' numbered';`);
  });

  await softStep('Apply Numbering Scheme dialog (engine + Bio applying it)', async () => {
    await page.evaluate(() => DG.Func.find({package: 'Bio'})
      .find((f) => f.name.replace(/\s/g, '').toLowerCase() === 'applynumberingscheme')!.apply({}));
    const dlg = page.locator('[name="dialog-Apply-Antibody-Numbering"]');
    await dlg.waitFor({timeout: 60_000});
    const t0 = Date.now();
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('AntibodyHC (aligned)'), null,
      {timeout: 1_800_000, polling: 500});
    await record(page, 'numbering dialog: OK to aligned column', Date.now() - t0);
  });

  fs.writeFileSync(`${OUT}/perf.json`, JSON.stringify({rows: ROWS, results}, null, 2));
  for (const e of stepErrors) console.log(`STEP ERROR ${e.step}: ${e.error}`);
});
