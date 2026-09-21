/* The steps only Bio can define: its readiness, and the readings of its own results — a molfile
   without the isotope flag that breaks standardization downstream (GROK-15176), a pairwise
   alignment's shape. Everything generic — the top menu, dialogs, columns, viewers — is the
   library's (`grok-bdd list-steps`). */
import {expect, Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {readResult} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** Bio initializes on its first call (RDKit, the monomer libraries, the sequence helper — about
 * eight seconds on a fresh page) and the platform holds every call of the package until then;
 * a service getter is that hold made visible, so a feature's first command is not the one that
 * pays for it. Free once the package is up — provided the helper stays in the page: returned to
 * Node it is the RDKit module's 16 MB heap serialized to base64, ten seconds per call. */
export const bioInitialized = Given('the Bio package is initialized', async (page: Page) => {
  await page.evaluate(async () => { await grok.functions.call('Bio:getSeqHelper', {}); });
}, {tier: 'api', description: 'awaits Bio\'s init through the platform (a call of a Bio service getter returns once initBio has run)'});

/** V3000 atom lines (`M  V30 <idx> <symbol> x y z ...`) carrying `MASS=1` on anything but H or D. */
export const noIsotopeFlag = Then('the result should not carry an isotope flag on a heavy atom', async (page: Page) => {
  const flagged: string[] = await readResult(page, `String(value ?? '').split('\\n').filter((l) => /^M  V30 \\d+ [A-Za-z]+ /.test(l) && / MASS=1\\b/.test(l) && !/^M  V30 \\d+ [HD] /.test(l))`);
  expect(flagged, 'V3000 heavy-atom lines with MASS=1 (GROK-15176)').toEqual([]);
}, {description: 'the last result read as a V3000 molfile: no heavy atom with MASS=1 — the flag that makes PubChem standardization reject it (GROK-15176)'});

/** A pairwise alignment result: its string fields are the aligned sequences, gap-padded to at
 * least the longer input. */
export const alignmentLength = Then('the result should be an alignment of at least {int} positions', async (page: Page, length: number) => {
  const lengths: number[] = await readResult(page, `Object.values(value ?? {}).filter((v) => typeof v === 'string').map((s) => s.length)`);
  expect(lengths.length, 'aligned sequences in the result').toBeGreaterThanOrEqual(2);
  expect(Math.min(...lengths), 'the shortest aligned sequence').toBeGreaterThanOrEqual(length);
});

/** The symbols of a HELM column: every `[Multi]` and every single letter inside the `{...}`
 * polymers, the way the column spells them. */
export const helmMonomersMatch = Then('the result should be exactly the monomers of {string} column', async (page: Page, column: string) => {
  const {result, parsed} = await page.evaluate((c) => {
    const col = grok.shell.t.col(c);
    if (!col)
      throw new Error(`no "${c}" column in the current table`);
    const symbols = new Set<string>();
    for (let i = 0; i < col.length; i++) {
      for (const block of String(col.get(i) ?? '').matchAll(/\{([^}]*)\}/g)) {
        for (const m of block[1].matchAll(/\[([^\]]+)\]|(?<![\w\[])([A-Za-z])(?![\w\]])/g))
          symbols.add(m[1] ?? m[2]);
      }
    }
    const value = (window as any).__bddLastResult?.value;
    return {result: Array.from(value ?? []).map(String).filter((s) => s !== '').sort(), parsed: [...symbols].sort()};
  }, column);
  expect(result, `the monomers returned against those in "${column}"`).toEqual(parsed);
});
