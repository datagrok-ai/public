import {test, expect, Page} from '@playwright/test';
import {baseUrl, specTestOptions, waitForMolecule, waitForChemMenu} from '../spec-login';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

test.use(specTestOptions);

const SPGI_FILE_URL = '/file/System.DemoFiles/SPGI.csv?browse=files';
const MOLBLOCK_MARKER = 'M  END';
const SAMPLE_ROWS = 25;

async function loginAndOpenSpgi(page: Page): Promise<void> {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test`, which derives the token from ~/.grok/config.yaml.');
  await page.goto(baseUrl + '/oauth/');
  await page.context().addCookies([{name: 'auth', value: token, domain: new URL(baseUrl).hostname, path: '/'}]);
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  await page.goto(baseUrl + SPGI_FILE_URL);
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout: 120_000});
}

test('Molecule column in exported CSV', async ({page}) => {
  test.setTimeout(300_000);

  await loginAndOpenSpgi(page);

  await waitForMolecule(page, 60_000);
  await waitForChemMenu(page);

  const result = await page.evaluate(async (sampleRows: number) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const df = grok.shell.tv.dataFrame;
    const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
    if (!molCol) throw new Error('No Molecule-semType column on SPGI');
    const extras = df.columns.toList()
      .filter((c: any) => c.semType !== 'Molecule')
      .slice(0, 2)
      .map((c: any) => c.name);
    const selected = [molCol.name, ...extras];

    const mask = DG.BitSet.create(df.rowCount, (i: number) => i < sampleRows);
    const small = df.clone(mask, selected);

    const smilesCsv = await small.toCsvEx({moleculesAsSmiles: true});
    const rawCsv = await small.toCsvEx({}); 
    return {
      smilesCsv, rawCsv, selected, mol: molCol.name,
      all: df.columns.names() as string[],
    };
  }, SAMPLE_ROWS);

  const lines = result.smilesCsv.split(/\r?\n/).filter((l) => l.length > 0);
  const header = parseCsvLine(lines[0]);

  expect(result.smilesCsv.length).toBeGreaterThan(0);
  expect(result.smilesCsv).not.toContain(MOLBLOCK_MARKER);

  expect(new Set(header)).toEqual(new Set(result.selected));
  for (const unselected of result.all.filter((n) => !result.selected.includes(n)))
    expect(header).not.toContain(unselected);

  const molIdx = header.indexOf(result.mol);
  expect(molIdx).toBeGreaterThanOrEqual(0);
  const firstRow = parseCsvLine(lines[1]);
  const firstSmiles = firstRow[molIdx];
  expect(firstSmiles.length).toBeGreaterThan(0);
  expect(firstSmiles).not.toContain(MOLBLOCK_MARKER);
  expect(firstSmiles).toMatch(/^[A-Za-z0-9@+\-\[\]\(\)=#%./\\:]+$/);

  expect(result.rawCsv).toContain(MOLBLOCK_MARKER);

  await test.step('Real CSV file is downloaded and verified on disk', async () => {
    const downloadPromise = page.waitForEvent('download', {timeout: 30_000});
    await page.evaluate(async (sampleRows: number) => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      const df = grok.shell.tv.dataFrame;
      const molCol = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const extras = df.columns.toList()
        .filter((c: any) => c.semType !== 'Molecule').slice(0, 2).map((c: any) => c.name);
      const small = df.clone(DG.BitSet.create(df.rowCount, (i: number) => i < sampleRows),
        [molCol.name, ...extras]);
      const csv = await small.toCsvEx({moleculesAsSmiles: true});
      DG.Utils.download((df.name || 'data') + '.csv', csv, 'text/csv');
    }, SAMPLE_ROWS);

    const download = await downloadPromise;
    expect(download.suggestedFilename()).toMatch(/\.csv$/i);

    const target = path.join(os.tmpdir(), `tt-spgi-export-${Date.now()}.csv`);
    await download.saveAs(target);
    try {
      const onDisk = fs.readFileSync(target, 'utf8');
      expect(onDisk.length).toBeGreaterThan(0);

      expect(onDisk).not.toContain(MOLBLOCK_MARKER);
      const dlHeader = parseCsvLine(onDisk.split(/\r?\n/)[0]);
      expect(new Set(dlHeader)).toEqual(new Set(result.selected));
    } finally {
      try { fs.unlinkSync(target); } catch (_) {  }
    }
  });
});

function parseCsvLine(line: string): string[] {
  const out: string[] = [];
  let cur = '';
  let inQ = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (inQ) {
      if (ch === '"') {
        if (line[i + 1] === '"') { cur += '"'; i++; }
        else inQ = false;
      } else cur += ch;
    } else if (ch === '"') inQ = true;
    else if (ch === ',') { out.push(cur); cur = ''; }
    else cur += ch;
  }
  out.push(cur);
  return out;
}
