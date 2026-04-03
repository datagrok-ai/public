import {test, expect, Page} from '@playwright/test';

const SPGI_FILE = 'System:AppData/Chem/tests/spgi-100.csv';
const PEPTIDES_FILE = 'System:DemoFiles/bio/peptides.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

async function login(page: Page) {
  await page.goto('/');
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
  });
  await page.waitForFunction(() => {
    const bar = document.querySelector('.grok-loader,.d4-loading-bar,.d4-progress-bar');
    return !bar || getComputedStyle(bar).display === 'none';
  }, {timeout: 60000}).catch(() => {});
  await page.waitForTimeout(3000);
}

async function openTable(page: Page, path: string, name: string) {
  await page.evaluate(async ({path, name}) => {
    grok.shell.closeAll();
    const df = await grok.data.files.openTable(path);
    df.name = name;
    grok.shell.addTableView(df);
  }, {path, name});
  await page.waitForTimeout(8000);
}

async function openFilterPanel(page: Page) {
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.waitForTimeout(3000);
}

async function getFilteredCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

test.describe('Chem and Bio: Filter Panel', () => {
  test.setTimeout(600_000);

  test('Chem and Bio filter scenario', async ({page}) => {
    stepErrors.length = 0;
    await login(page);

    // #### Section 1: Refresh with chem filter

    await softStep('1.1 Open SPGI dataset', async () => {
      await openTable(page, SPGI_FILE, 'SPGI');
      const totalRows = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(totalRows).toBeGreaterThan(0);
    });

    await softStep('1.2 Open the Filter Panel', async () => {
      await openFilterPanel(page);
      const hasStructure = await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        return fg.filters.some((f: any) => f.column && f.column.name === 'Structure');
      });
      expect(hasStructure).toBe(true);
    });

    await softStep('1.3 Draw CCC(N(C)C)=O in Structure filter', async () => {
      // Click Sketch link to open sketcher
      const sketchLink = await page.evaluate(() => {
        const links = document.querySelectorAll('.sketch-link, .d4-link-action');
        for (const l of links) {
          if (l.textContent!.trim() === 'Sketch') {
            const card = l.closest('.d4-filter');
            if (card) {
              const header = card.querySelector('.d4-filter-column-name');
              if (header && header.textContent!.trim() === 'Structure') {
                (l as HTMLElement).click();
                return true;
              }
            }
          }
        }
        // Fallback: click the first Sketch link (Structure is first chem filter)
        if (links.length > 0) {
          (links[0] as HTMLElement).click();
          return true;
        }
        return false;
      });
      await page.waitForTimeout(2000);

      // Type SMILES in the sketcher dialog
      const smilesInput = page.getByPlaceholder('SMILES, MOLBLOCK, Inchi, ChEMBL id, etc');
      if (await smilesInput.isVisible({timeout: 5000}).catch(() => false)) {
        await smilesInput.fill('CCC(N(C)C)=O');
        await page.keyboard.press('Enter');
        await page.waitForTimeout(1000);
        await page.locator('button:has-text("OK")').click();
      }

      // Wait for substructure search to complete
      await page.waitForTimeout(10000);
      const filtered = await getFilteredCount(page);
      const total = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(filtered).toBeLessThan(total);
      expect(filtered).toBeGreaterThan(0);
    });

    await softStep('1.4 Switch search types and verify row counts', async () => {
      const results = await page.evaluate(async () => {
        const fg = grok.shell.tv.getFiltersGroup();
        let structFilter: any = null;
        for (const f of fg.filters)
          if (f.column && f.column.name === 'Structure') { structFilter = f; break; }
        if (!structFilter) return null;

        const types = ['Included in', 'Exact', 'Similar', 'Not contains', 'Not included in', 'Contains'];
        const results: Record<string, number> = {};
        for (const t of types) {
          structFilter.searchTypeInput.value = t;
          structFilter.searchTypeInput.fireChanged();
          for (let i = 0; i < 30; i++) {
            await new Promise(r => setTimeout(r, 1000));
            if (!structFilter.searchNotCompleted) break;
          }
          results[t] = grok.shell.tv.dataFrame.filter.trueCount;
        }
        return results;
      });

      expect(results).not.toBeNull();
      if (results) {
        // Included in, Exact, Similar should be 0 (molecule not in dataset)
        expect(results['Included in']).toBe(0);
        expect(results['Exact']).toBe(0);
        expect(results['Similar']).toBe(0);
        // Not contains + Contains should equal total rows
        const total = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
        expect(results['Not contains']).toBeGreaterThan(0);
        expect(results['Not included in']).toBe(total);
        expect(results['Contains']).toBeGreaterThan(0);
      }
    });

    await softStep('1.5 Close All', async () => {
      await page.evaluate(() => grok.shell.closeAll());
      await page.waitForTimeout(1000);
    });

    // #### Section 2: Bio

    await softStep('2.1 Open peptides.csv', async () => {
      await openTable(page, PEPTIDES_FILE, 'peptides');
      const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
      expect(rowCount).toBe(647);
    });

    await softStep('2.2 Wait for AlignedSequence to render as Macromolecule', async () => {
      await page.waitForTimeout(5000);
      const semType = await page.evaluate(() =>
        grok.shell.tv.dataFrame.col('AlignedSequence').semType
      );
      expect(semType).toBe('Macromolecule');
    });

    await softStep('2.3 Open the Filter Panel', async () => {
      await openFilterPanel(page);
      const hasAlignedSeq = await page.evaluate(() => {
        const fg = grok.shell.tv.getFiltersGroup();
        return fg.filters.some((f: any) => f.column && f.column.name === 'AlignedSequence');
      });
      expect(hasAlignedSeq).toBe(true);
    });

    await softStep('2.4 Enter T-T-Y-K-N-Y-V in Substructure filter', async () => {
      const filtered = await page.evaluate(async () => {
        const fg = grok.shell.tv.getFiltersGroup();
        let bioFilter: any = null;
        for (const f of fg.filters)
          if (f.column && f.column.name === 'AlignedSequence') { bioFilter = f; break; }
        if (!bioFilter || !bioFilter.bioFilter) return -1;

        bioFilter.bioFilter.substructureInput.value = 'T-T-Y-K-N-Y-V';
        bioFilter.bioFilter.substructureInput.fireChanged();
        await new Promise(r => setTimeout(r, 5000));
        return grok.shell.tv.dataFrame.filter.trueCount;
      });
      expect(filtered).toBeGreaterThan(0);
      expect(filtered).toBeLessThan(647);
    });

    await softStep('2.5 Close All', async () => {
      await page.evaluate(() => grok.shell.closeAll());
      await page.waitForTimeout(1000);
    });

    if (stepErrors.length > 0) {
      const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
      throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
    }
  });
});
