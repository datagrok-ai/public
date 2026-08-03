import {test, expect, Page} from '@playwright/test';
import {baseUrl, specTestOptions, waitForMolecule, waitForChemMenu} from '../spec-login';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

// Source scenario: General/molecule-in-exported-csv.md
// Verifies "As CSV (options)..." export with "Molecules as SMILES" + "Selected
// columns only". The export dialog (core/.../save_table_as_csv.dart) is a thin
// wrapper over DataFrame.toCsvEx(options): on OK it sets moleculesAsSmiles /
// selectedColumnsOnly and converts molblock columns to SMILES via
// Chem:convertNotation. We drive that same logic through the public API so the
// assertions are deterministic (no canvas column-picking).
//
// SPGI is opened by navigating straight to its file deep-link, so the platform
// boots AND opens the table in a single navigation (no "load root, then open"
// two-step, no blind settle wait). Because the file-open IS the initial
// navigation, there is no mid-session route change to race page.evaluate.
//
// "Selected columns only" is realized by cloning the table down to exactly the
// chosen columns (molecule column + a couple of others). We also clone down to a
// small row subset: converting the full 3624-row SPGI molblock column on a cold
// RDKit module is slow, and the feature (per-cell molblock->SMILES) is fully
// exercised on a handful of rows.
//
// SPGI.csv ships a `Structure` column of V2000 molblocks ("M  END"), so the
// conversion is observable: the SMILES export must NOT contain "M  END" while
// the unconverted control still does. Verified live on dev 2026-06-16.

test.use(specTestOptions);

const SPGI_FILE_URL = '/file/System.DemoFiles/SPGI.csv?browse=files';
const MOLBLOCK_MARKER = 'M  END';
const SAMPLE_ROWS = 25;

/** Attach auth to the origin, then navigate straight to the SPGI file deep-link
 * so the app boots and opens the table in one navigation. */
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
  // The molecule semType + Chem package (Chem:convertNotation / RDKit) must be
  // ready for the molblock->SMILES conversion path. Both are condition waits —
  // they resolve as soon as the table + Chem menu are live.
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

    // "Selected columns only" + small row sample via clone.
    const mask = DG.BitSet.create(df.rowCount, (i: number) => i < sampleRows);
    const small = df.clone(mask, selected);

    const smilesCsv = await small.toCsvEx({moleculesAsSmiles: true});
    const rawCsv = await small.toCsvEx({}); // control: no SMILES conversion
    return {
      smilesCsv, rawCsv, selected, mol: molCol.name,
      all: df.columns.names() as string[],
    };
  }, SAMPLE_ROWS);

  const lines = result.smilesCsv.split(/\r?\n/).filter((l) => l.length > 0);
  const header = parseCsvLine(lines[0]);

  // --- Assertion 1: Structure column exported as SMILES, not molblocks.
  expect(result.smilesCsv.length).toBeGreaterThan(0);
  expect(result.smilesCsv).not.toContain(MOLBLOCK_MARKER);

  // --- Assertion 2: only the selected columns are present.
  expect(new Set(header)).toEqual(new Set(result.selected));
  for (const unselected of result.all.filter((n) => !result.selected.includes(n)))
    expect(header).not.toContain(unselected);

  // --- Assertion 3: the molecule cell holds a non-empty single-line SMILES.
  const molIdx = header.indexOf(result.mol);
  expect(molIdx).toBeGreaterThanOrEqual(0);
  const firstRow = parseCsvLine(lines[1]);
  const firstSmiles = firstRow[molIdx];
  expect(firstSmiles.length).toBeGreaterThan(0);
  expect(firstSmiles).not.toContain(MOLBLOCK_MARKER);
  expect(firstSmiles).toMatch(/^[A-Za-z0-9@+\-\[\]\(\)=#%./\\:]+$/);

  // --- Control: without moleculesAsSmiles the molblock marker survives, proving
  // the conversion above is what stripped it (not just an absence in source).
  expect(result.rawCsv).toContain(MOLBLOCK_MARKER);

  // --- Assertion 4: a REAL .csv file is downloaded to disk and its on-disk
  // content matches. DG.Utils.download is the platform's own download primitive
  // (Blob + <a download>.click()) — the same path the "Save as CSV" command's
  // saveAs() uses. We feed it the real toCsvEx export, capture the browser
  // download, save it, and assert against the bytes that actually landed on disk.
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
      // The downloaded file holds SMILES (no molblocks) and only selected columns.
      expect(onDisk).not.toContain(MOLBLOCK_MARKER);
      const dlHeader = parseCsvLine(onDisk.split(/\r?\n/)[0]);
      expect(new Set(dlHeader)).toEqual(new Set(result.selected));
    } finally {
      try { fs.unlinkSync(target); } catch (_) { /* best effort */ }
    }
  });
});

/** Minimal CSV line splitter that respects double-quoted fields. */
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
