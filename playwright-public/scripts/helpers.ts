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
  // Wait for #grok-preloader to detach. On a cold CI Datlas the Browse tree
  // becomes visible BEFORE the preloader is dismissed; while present, it covers
  // rootDiv and intercepts clicks (incl. on dropdown menu items like "Pyodide
  // Script..."), surfacing as "subtree intercepts pointer events" timeouts.
  await page.waitForFunction(
    () => document.querySelector('#grok-preloader, .grok-preloader') == null,
    undefined, { timeout: 30_000 },
  ).catch(() => { /* tolerate a lingering preloader — neutralised by addStyleTag below */ });
  // Also neutralise pointer events on the preloader — it can re-appear after
  // dialog-opening clicks while the platform fetches data, and intercept the
  // next interaction. Platform JS still drives it; only the click-intercept
  // is disabled for the page lifetime. The same `display: none` for
  // `.d4-tooltip` mirrors connections/queries helpers — hover-tooltips
  // (e.g. the "After you save the project ..." hint that pops up next to
  // the Save-project dialog's OK button) can also overlap actionability.
  await page.addStyleTag({ content: `
    #grok-preloader, .grok-preloader { pointer-events: none !important; }
    .d4-tooltip { display: none !important; }
  ` });

  // Expand Platform node (only if collapsed)
  const platformExpander = page.locator('[name="tree-expander-Platform"]');
  await expect(platformExpander).toBeVisible({ timeout: 15_000 });
  if (!(await platformExpander.evaluate(el => el.classList.contains('d4-tree-view-tri-expanded')))) {
    await platformExpander.click();
    await expect(platformExpander).toHaveClass(/d4-tree-view-tri-expanded/, { timeout: 5_000 });
  }

  // Expand Functions node (only if collapsed)
  const functionsExpander = page.locator('[name="tree-expander-Platform---Functions"]');
  await expect(functionsExpander).toBeVisible({ timeout: 5_000 });
  if (!(await functionsExpander.evaluate(el => el.classList.contains('d4-tree-view-tri-expanded')))) {
    await functionsExpander.click();
    await expect(functionsExpander).toHaveClass(/d4-tree-view-tri-expanded/, { timeout: 5_000 });
  }

  // Click Scripts label
  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  await expect(scriptsLabel).toBeVisible({ timeout: 5_000 });
  await scriptsLabel.click();
  await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 15_000 });
  // Wait for the gallery to actually render at least one card — ensures data
  // is loaded before tests start searching. Without this, searching immediately
  // after the search bar appears can filter against an empty dataset.
  await page.locator('.grok-gallery-grid-item').first().waitFor({ state: 'visible', timeout: 15_000 }).catch(() => { /* fresh CI stack has no scripts yet — gallery legitimately empty; tests create their own below */ });
}

// Close all open views and return to the Scripts gallery without a full page
// reload. Fast reset between tests in shared-context suites.
export async function resetToScripts(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    if (g?.shell?.closeAll) g.shell.closeAll();
    document.querySelectorAll('.d4-dialog').forEach((d: any) => { try { d.remove(); } catch (_) {} });
    document.querySelectorAll('.d4-toast, .d4-balloon, .d4-menu').forEach((e) => e.remove());
  });
  const scriptsLabel = page.locator('.d4-tree-view-item-label', { hasText: /^Scripts$/i }).first();
  if (await scriptsLabel.isVisible({ timeout: 1_000 }).catch(() => false)) {
    await scriptsLabel.click();
    await expect(page.locator('.grok-gallery-search-bar')).toBeVisible({ timeout: 10_000 });
  } else {
    await openScriptsBrowser(page);
  }
  // Gallery cards are populated asynchronously after the search bar appears.
  // Without this wait, an immediate searchScript() call filters against an
  // empty dataset (debounced search caches "no results") and times out.
  await page.locator('.grok-gallery-grid-item').first().waitFor({ state: 'visible', timeout: 15_000 }).catch(() => { /* fresh CI stack has no scripts yet — gallery legitimately empty; tests create their own below */ });
}

// Idempotent: create the testRscript only if it does not exist.
// Cheaper than delete+create in beforeEach when the script can persist across tests.
export async function apiEnsureScript(page: Page) {
  await page.evaluate(async ([name, content]) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    const cap = name.charAt(0).toUpperCase() + name.slice(1);
    for (const q of [name, cap]) {
      const list = await grok.dapi.scripts.filter(`name = "${q}"`).list();
      if (list.length > 0) return;
    }
    const script = DG.Script.create(content as string);
    await grok.dapi.scripts.save(script);
  }, [SCRIPT_NAME, R_SCRIPT_CONTENT] as [string, string]);
}

// Close all open views and transient UI via shell.closeAll — returns the
// application to its initial condition between tests without a full page
// reload. Also clears leftover dialogs, toasts, and context menus.
export async function resetShell(page: Page) {
  await page.evaluate(() => {
    const g = (window as any).grok;
    if (g?.shell?.closeAll) g.shell.closeAll();
    document.querySelectorAll('.d4-dialog').forEach((d: any) => { try { d.remove(); } catch (_) {} });
    document.querySelectorAll('.d4-toast, .d4-balloon, .d4-menu').forEach((e) => e.remove());
  });
  await page.waitForTimeout(300);
}

// Search the scripts gallery for a script name and return the card locator.
// Waits for the card to render after the debounced search filter runs.
export async function searchScript(page: Page, name: string) {
  const searchInput = page.locator('input[placeholder="Search scripts by name or by #tags"]');
  await expect(searchInput).toBeVisible();
  await searchInput.fill('');
  await searchInput.fill(name);
  await searchInput.press('Enter');
  const card = getScriptCard(page, name);
  await expect(card).toBeVisible({ timeout: 20_000 });
  return card;
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
