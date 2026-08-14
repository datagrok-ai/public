/* ---
sub_features_covered: [bio.analyze.activity-cliffs.top-menu, bio.analyze.compare-sequences.top-menu, bio.analyze.composition.top-menu, bio.analyze.msa.top-menu, bio.analyze.sequence-space.top-menu, bio.annotate.manage-annotations.top-menu, bio.annotate.numbering-scheme.top-menu, bio.annotate.scan-liabilities.top-menu, bio.calculate.get-region.top-menu, bio.calculate.identity.top-menu, bio.calculate.similarity.top-menu, bio.manage.libraries-view.top-menu, bio.manage.match-with-library.top-menu, bio.manage.monomers-view.top-menu, bio.search.diversity.top-menu, bio.search.similarity.top-menu, bio.search.subsequence.top-menu, bio.top-menu.registration, bio.transform.convert-notation.top-menu, bio.transform.molecules-to-helm.top-menu, bio.transform.split-to-monomers.top-menu, bio.transform.to-atomic-level.top-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {openBioGroup, openBioMenu, dialog, cancelDialog, loadBioTable} from './helpers';

test.use(specTestOptions);

const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const TABLE = 'antibodies';
const ROWS = 40;

/** What a Bio top-menu leaf is expected to produce. Exactly one of dialog/viewer/view is set —
 * the dialog `name` attribute is the function's friendly name, which is not always the menu label. */
interface Leaf {
  group: string;
  leaf: string;
  dialog?: string;
  viewer?: string;
  view?: string;
}

const BIO_LEAVES: Leaf[] = [
  {group: 'Transform', leaf: 'Molecules-to-HELM...', dialog: 'Molecules-to-HELM'},
  {group: 'Transform', leaf: 'To-Atomic-Level...', dialog: 'To-Atomic-Level'},
  {group: 'Transform', leaf: 'Convert-Sequence-Notation...', dialog: 'Convert-Sequence-Notation'},
  {group: 'Transform', leaf: 'Split-to-Monomers...', dialog: 'Split-to-Monomers'},
  {group: 'Analyze', leaf: 'Activity-Cliffs...', dialog: 'Sequence-Activity-Cliffs'},
  {group: 'Analyze', leaf: 'Sequence-Space...', dialog: 'Sequence-Space'},
  {group: 'Analyze', leaf: 'MSA...', dialog: 'MSA'},
  {group: 'Analyze', leaf: 'Compare-sequences...', dialog: 'Compare-Sequences'},
  {group: 'Analyze', leaf: 'Composition', dialog: 'Composition-Analysis'},
  {group: 'Calculate', leaf: 'Extract-Region...', dialog: 'Get-Sequence-Region'},
  {group: 'Calculate', leaf: 'Identity...', dialog: 'Identity'},
  {group: 'Calculate', leaf: 'Similarity...', dialog: 'Similarity'},
  {group: 'Annotate', leaf: 'Apply-Numbering-Scheme...', dialog: 'Apply-Antibody-Numbering'},
  {group: 'Annotate', leaf: 'Scan-Liabilities...', dialog: 'Scan-Sequence-Liabilities'},
  {group: 'Annotate', leaf: 'Manage-Annotations...', dialog: 'Manage-Annotations'},
  {group: 'Manage', leaf: 'Match-with-Monomer-Library...', dialog: 'matchWithMonomerLibrary'},
  {group: 'Manage', leaf: 'Monomer-Libraries', view: 'Manage Monomer Libraries'},
  {group: 'Manage', leaf: 'Monomers', view: 'Manage Monomers'},
  {group: 'Search', leaf: 'Similarity-Search', viewer: 'Sequence Similarity Search'},
  {group: 'Search', leaf: 'Diversity-Search', viewer: 'Sequence Diversity Search'},
  {group: 'Search', leaf: 'Subsequence-Search-...', dialog: 'Substructure-Search'},
];

/** Leaves other packages contribute to the Bio menu. They are smoke-checked only when their
 * package happens to be installed — Bio's suite must not demand SequenceTranslator or Peptides.
 * Bio | Folding is deliberately absent: EsmFold/Boltz submit real inference jobs. */
const FOREIGN_LEAVES: Leaf[] = [
  {group: 'Analyze', leaf: 'Hierarchical-Clustering...', dialog: 'Hierarchical-Clustering'},
  {group: 'Analyze', leaf: 'SAR...', dialog: 'Analyze-Peptides'},
  {group: 'Transform', leaf: 'Fetch-PDB-Sequences...', dialog: 'Fetch-PDB-Sequences'},
  {group: 'PolyTool', leaf: 'Enumerate-HELM...', dialog: 'PolyTool-Helm-Enumeration'},
  {group: 'PolyTool', leaf: 'Combine-Sequences...', dialog: 'Combine-Sequences'},
];

const label = (l: Leaf): string => `Bio | ${l.group} | ${l.leaf.replace(/-/g, ' ').trim()}`;

/** Returns to the table view and waits until the Bio menu is usable again. */
async function backToTable(page: Page): Promise<void> {
  await page.evaluate((n) => {
    const tv = Array.from(grok.shell.tableViews).find((v: any) => v.name === n);
    if (tv && grok.shell.v !== tv) grok.shell.v = tv as any;
  }, TABLE);
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 30_000});
}

test('Bio top menu — every Bio-owned leaf is registered under its group', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await loadBioTable(page, ANTIBODIES, ROWS, TABLE);

  for (const group of ['Transform', 'Analyze', 'Calculate', 'Annotate', 'Manage', 'Search']) {
    await softStep(`Bio | ${group} lists its leaves`, async () => {
      await openBioGroup(page, group);
      for (const l of BIO_LEAVES.filter((x) => x.group === group))
        await expect(page.locator(`[name="div-Bio---${group}---${l.leaf}"]`)).toBeVisible({timeout: 15_000});
    });
  }
  await page.keyboard.press('Escape');
  finishSpec();
});

test('Bio top menu — every leaf opens its dialog, viewer or view and closes cleanly', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await loadBioTable(page, ANTIBODIES, ROWS, TABLE);

  for (const l of BIO_LEAVES) {
    await softStep(`${label(l)} opens and closes`, async () => {
      await backToTable(page);
      await openBioMenu(page, l.group, l.leaf);

      if (l.dialog) {
        const dlg = dialog(page, l.dialog);
        await dlg.waitFor({timeout: 90_000});
        // A dialog with no OK cannot be run at all — that is a broken entry point, not a slow one.
        expect(await dlg.locator('[name="button-OK"]').count()).toBe(1);
        await cancelDialog(page, l.dialog);
      }
      else if (l.viewer) {
        await page.waitForFunction((t) => Array.from(((grok.shell.tv as any)?.viewers ?? []) as any[])
          .some((v: any) => v.type === t), l.viewer, {timeout: 120_000});
        await page.evaluate((t) => {
          for (const v of Array.from(((grok.shell.tv as any)?.viewers ?? []) as any[]))
            if (v.type === t) v.close();
        }, l.viewer);
      }
      else {
        await page.waitForFunction((n) => grok.shell.v?.name === n, l.view, {timeout: 120_000});
        await page.evaluate((n) => { if (grok.shell.v?.name === n) grok.shell.v.close(); }, l.view);
      }

      expect(await page.locator('.d4-dialog').count()).toBe(0);
    });
  }

  await backToTable(page);
  finishSpec();
});

test('Bio top menu — leaves contributed by other packages open when those packages are installed',
  async ({page}) => {
    test.setTimeout(600_000);
    stepErrors.length = 0;
    await loginToDatagrok(page);
    await loadBioTable(page, ANTIBODIES, ROWS, TABLE);

    const skipped: string[] = [];
    for (const l of FOREIGN_LEAVES) {
      await softStep(`${label(l)} (foreign) opens and closes`, async () => {
        await backToTable(page);
        await page.locator('[name="div-Bio"]').click();
        const groupLoc = page.locator(`[name="div-Bio---${l.group}"]`);
        if (await groupLoc.count() === 0) {
          skipped.push(label(l));
          await page.keyboard.press('Escape');
          return;
        }
        await openBioGroup(page, l.group);
        const leafLoc = page.locator(`[name="div-Bio---${l.group}---${l.leaf}"]`);
        if (await leafLoc.count() === 0) {
          skipped.push(label(l));
          await page.keyboard.press('Escape');
          return;
        }
        await leafLoc.click();
        await dialog(page, l.dialog!).waitFor({timeout: 90_000});
        await cancelDialog(page, l.dialog!);
      });
    }
    if (skipped.length > 0)
      console.log(`[note] not installed on this stand, so not exercised: ${skipped.join(', ')}`);
    finishSpec();
  });
