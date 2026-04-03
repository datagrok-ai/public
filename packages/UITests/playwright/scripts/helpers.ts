import { Page, expect } from '@playwright/test';

const BASE = process.env.DATAGROK_URL!;

export const SCRIPT_NAME = 'testRscript';

// Full R script content used across create/edit/run tests
export const R_SCRIPT_CONTENT = `#name: ${SCRIPT_NAME}
#language: r
#sample: cars.csv
#input: dataframe table [Data table]
#output: int count [Number of cells in table]
#output: string newParam

count <- nrow(table) * ncol(table)
newParam <- "test"`;

// Navigate to Scripts via Browse > Platform > Functions > Scripts
export async function openScriptsBrowser(page: Page) {
  await page.goto(`${BASE}/browse`);

  // Expand Platform node (only if collapsed)
  const platformExpander = page.locator('[name="tree-expander-Platform"]');
  await expect(platformExpander).toBeVisible({ timeout: 15_000 });
  if (!(await platformExpander.evaluate(el => el.classList.contains('d4-tree-view-tri-expanded')))) {
    await platformExpander.click();
    await page.waitForTimeout(1000);
  }

  // Expand Functions node (only if collapsed)
  const functionsExpander = page.locator('[name="tree-expander-Platform---Functions"]');
  await expect(functionsExpander).toBeVisible({ timeout: 5_000 });
  if (!(await functionsExpander.evaluate(el => el.classList.contains('d4-tree-view-tri-expanded')))) {
    await functionsExpander.click();
    await page.waitForTimeout(1000);
  }

  // Click Scripts label
  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  await expect(scriptsLabel).toBeVisible({ timeout: 5_000 });
  await scriptsLabel.click();
  await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 15_000 });
}


// Find a script card by name in the gallery
export function getScriptCard(page: Page, name: string) {
  return page.locator('.grok-gallery-grid-item', {
    has: page.locator('.grok-gallery-grid-item-title', { hasText: name }),
  }).first();
}

// Right-click a script card to open its context menu.
// Uses 'Delete' as a sentinel because it always appears in the context menu for
// user-owned scripts and cannot be confused with the hidden ribbon menu items.
export async function rightClickScript(page: Page, name: string) {
  const card = getScriptCard(page, name);
  await expect(card).toBeVisible({ timeout: 10_000 });
  await card.dispatchEvent('contextmenu', { button: 2, bubbles: true });
  await expect(page.locator('.d4-menu-item-label', { hasText: 'Delete' }).first()).toBeVisible({ timeout: 5_000 });
}

// Click a context menu item by its visible label
export async function clickMenuItem(page: Page, label: string) {
  await page.locator('.d4-menu-item-label', { hasText: label }).first().click();
  await page.waitForTimeout(600);
}

// Set CodeMirror editor content in the script editor.
// Uses keyboard input (Ctrl+A, then type) instead of cm.setValue() because
// programmatic setValue() does not trigger Datagrok's dirty-state detection.
export async function setScriptContent(page: Page, content: string) {
  await page.locator('.CodeMirror').click();
  await page.keyboard.press('ControlOrMeta+a');
  await page.keyboard.type(content, { delay: 0 });
  await page.waitForTimeout(400);
}

// Delete a script by name via API (for test setup/teardown).
// Datagrok capitalizes the name (testRscript → TestRscript), so search both variants.
export async function apiDeleteScript(page: Page, name: string) {
  await page.evaluate(async (n) => {
    try {
      const grok = (window as any).grok;
      const cap = n.charAt(0).toUpperCase() + n.slice(1);
      for (const q of [n, cap]) {
        const list = await grok.dapi.scripts.filter(`name = "${q}"`).list();
        for (const s of list) await grok.dapi.scripts.delete(s);
      }
    } catch (_) {
      // Script may not exist yet — ignore
    }
  }, name);
}

// Create the testRscript via API (prerequisite for edit/run/browse/delete tests).
// Uses DG.Script.create() which parses annotation headers from the code string.
// Returns null explicitly to prevent Playwright serialization errors with complex Dart objects.
export async function apiCreateScript(page: Page) {
  const created = await page.evaluate(async ([name, content]) => {
    const DG = (window as any).DG;
    const grok = (window as any).grok;
    const script = DG.Script.create(content as string);
    const saved = await grok.dapi.scripts.save(script);
    return saved?.name || saved?.friendlyName || 'unknown';
  }, [SCRIPT_NAME, R_SCRIPT_CONTENT] as [string, string]);
  await page.waitForTimeout(1000);
}

// Load cars.csv demo table (setup utility — uses API to load data reliably).
// This is a prerequisite for tests, not a tested UI action.
// Opens cars.csv via the Datagrok file API. Sets table name to "cars" for console compatibility.
export async function loadCarsDemoTable(page: Page) {
  await page.evaluate(() => {
    void (window as any).grok.dapi.files.readCsv('System:DemoFiles/cars.csv').then(
      (t: any) => { t.name = 'cars'; void (window as any).grok.shell.addTableView(t); },
      () => {
        void (window as any).grok.data.getDemoTable('cars.csv').then(
          (t2: any) => { t2.name = 'cars'; void (window as any).grok.shell.addTableView(t2); },
          () => {},
        );
      },
    );
  });
  await expect(page.locator('.d4-grid')).toBeVisible({ timeout: 20_000 });
}

// Close all open views by clicking the X button on each tab
export async function closeAllViews(page: Page) {
  const closeBtns = page.locator('.tab-handle-close');
  let count = await closeBtns.count();
  while (count > 0) {
    await closeBtns.first().click();
    await page.waitForTimeout(300);
    count = await closeBtns.count();
  }
  await page.waitForTimeout(300);
}

// Expand an accordion pane in the context panel by header text
export async function expandAccordionPane(page: Page, title: string) {
  const header = page.locator('.d4-accordion-pane-header', { hasText: title }).first();
  await expect(header).toBeVisible({ timeout: 8_000 });
  await header.click();
  await page.waitForTimeout(500);
}
