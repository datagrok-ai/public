import { test as base, expect, type Page } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

/**
 * Shared, already-authenticated page reused across every test in the worker.
 * `global-setup` logs in and writes `e2e/.auth.json`; we boot the Datagrok SPA
 * exactly once here, then tests load data and reset state through the JS API
 * (`openDemoCsv` / `resetShell` below) with no further full-page navigation —
 * this removes the two ~3s reloads that previously dominated each test.
 */
export const test = base.extend<{}, { sharedPage: Page }>({
  sharedPage: [async ({ browser }, use) => {
    const context = await browser.newContext({
      storageState: 'e2e/.auth.json',
      viewport: { width: 1920, height: 1080 },
    });
    const page = await context.newPage();
    await page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: 120_000 });
    await page.waitForFunction(() => document.querySelector('.grok-preloader') == null,
      undefined, { timeout: 120_000 });
    await page.locator('.d4-ribbon').first().waitFor({ timeout: 60_000 });
    await use(page);
    await context.close();
  }, { scope: 'worker' }],
  // Override the built-in test-scoped `page` so existing `async ({ page }) => …`
  // bodies transparently receive the shared worker page.
  page: async ({ sharedPage }, use) => {
    await use(sharedPage);
  },
});

export { expect };

/** `[name="input-host-..."]` wrapper around any DG input. */
export const inputHost = (safeName: string): string => `[name="input-host-${safeName}"]`;

/** Editable element inside an input host. */
export const inputEditor = (safeName: string): string => `${inputHost(safeName)} .ui-input-editor`;

/**
 * Open a CSV from `System:DemoFiles/` into a fresh TableView via the JS API
 * (`readCsv` + `addTableView`) — no page navigation. With the shared worker page
 * the SPA is already booted, so this is the ~0.4s fast path instead of a ~3s
 * Files-browser navigation + dblclick. The ML menu attaches to the new
 * TableView exactly as it does for a UI-opened file.
 *
 * `fileName` is the bare CSV name (e.g. `cars.csv`); the folder is always `System:DemoFiles`.
 */
export async function openDemoCsv(page: Page, fileName: string, uiTimeoutMs = 25_000): Promise<void> {
  await page.evaluate(async (name) => {
    const g = (window as unknown as { grok: any }).grok;
    const df = await g.dapi.files.readCsv(`System:DemoFiles/${name}`);
    g.shell.addTableView(df);
  }, fileName);
  await waitForCurrentTableView(page, uiTimeoutMs);
}

/** Resolves true once `grok.shell.tv?.dataFrame` is populated. */
export async function waitForCurrentTableView(page: Page, timeoutMs: number): Promise<boolean> {
  try {
    await page.waitForFunction(() => {
      const g = (window as unknown as { grok?: any }).grok;
      return !!g?.shell?.tv?.dataFrame?.columns?.length;
    }, undefined, { timeout: timeoutMs });
    return true;
  } catch {
    return false;
  }
}

/**
 * Click a Top-Menu leaf by its `name=` attribute. Top-menu items are tagged
 * `[name="div-A---B---C..."]` where the path uses `---` separators and trailing dots
 * are preserved from the original label (e.g. `PCA...` -> `div-ML---Analyze---PCA...`).
 *
 * If the leaf isn't rendered yet, walk the ancestor chain (`div-A`, `div-A---B`, ...)
 * priming each via click + hover events — nested submenus (e.g. `ML > Models > Train Model`)
 * require hovering the parent before the leaf is mounted. Direct DOM events are used
 * rather than `page.click` because the platform's menu handlers listen at the DOM level
 * and pointer-actionability checks fail intermittently on hover-revealed items.
 */
export async function clickTopMenuLeaf(page: Page, nameAttr: string): Promise<void> {
  const prefix = nameAttr.startsWith('div-') ? 'div-' : '';
  const path = nameAttr.replace(/^div-/, '').split('---');

  for (let depth = 1; depth < path.length; depth++) {
    const parentSel = `[name="${prefix}${path.slice(0, depth).join('---')}"]`;
    await page.evaluate((sel) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el) return;
      el.click();
      el.dispatchEvent(new MouseEvent('mouseenter', { bubbles: true }));
      el.dispatchEvent(new MouseEvent('mousemove', { bubbles: true }));
    }, parentSel);
    await page.waitForTimeout(350);
  }

  const ok = await page.evaluate((sel) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el) return false;
    el.click();
    return true;
  }, `[name="${nameAttr}"]`);
  if (!ok) throw new Error(`Top-menu leaf not found: ${nameAttr}`);
  await page.waitForTimeout(400);
}

/**
 * Wait for the next `.d4-dialog` whose title contains `title` (case-insensitive, substring).
 * Substring match is intentional — different ML dialogs render the title with mixed casing
 * (e.g. "Multivariate analysis" vs the menu's "Multivariate Analysis...") and we only need
 * to know that *some* matching dialog is mounted before we start filling fields.
 */
export async function waitForDialog(page: Page, title: string, timeoutMs = 15_000): Promise<void> {
  const re = new RegExp(escapeRegex(title), 'i');
  await page.locator('.d4-dialog .d4-dialog-title', { hasText: re }).first().waitFor({ timeout: timeoutMs });
}

/** Locator for the topmost open dialog. */
export function topDialog(page: Page) {
  return page.locator('.d4-dialog').last();
}

/**
 * Click the dialog's primary action. The platform uses inconsistent casing across ML
 * dialogs — ANOVA exposes `button-Run`, MVA/PLS expose `button-RUN`, PCA exposes
 * `button-OK`. Try all three so callers don't need to know which one their dialog uses.
 */
export async function clickDialogPrimary(page: Page, candidates: string[] = ['OK', 'Run', 'RUN']): Promise<void> {
  for (const c of candidates) {
    const btn = topDialog(page).locator(`[name="button-${c}"]`).first();
    if (await btn.count() > 0 && await btn.isEnabled().catch(() => false)) {
      await btn.click();
      return;
    }
  }
  throw new Error(`No enabled primary action found among: ${candidates.join(', ')}`);
}

/**
 * Open the column-picker sub-dialog for a multi-column input host and click `All` then OK.
 * Used by PCA/PLS scenarios that say "select all available columns" — the All link selects
 * every numeric column, which is the only target the picker offers for analytic inputs.
 */
export async function selectAllColumnsInPicker(page: Page, hostName: string): Promise<void> {
  await page.locator(`${inputHost(hostName)} .ui-input-editor`).first().click();
  // The picker is a fresh `.d4-dialog` layered on top of the analysis dialog.
  await page.locator('.d4-dialog [name="label-All"]').last().waitFor({ timeout: 10_000 });
  await page.locator('.d4-dialog [name="label-All"]').last().click();
  await page.locator('.d4-dialog [name="button-OK"]').last().click();
  await page.waitForTimeout(500);
}

/**
 * Replace the value of a multi-column input in the currently-open analysis dialog by reaching
 * through the dialog API: `DG.Dialog.getOpenDialogs()[-1].input(caption).value = cols`.
 *
 * Used as the deterministic fallback when the canvas-rendered "Select columns..." sub-dialog
 * cannot toggle individual rows via DOM events (see pls-run.md and linear-regression-run.md
 * retrospectives). The walkthrough — menu navigation, Components input, RUN — stays in the UI;
 * only this single field-set bypasses the unreachable canvas picker.
 */
export async function setDialogColumnListInput(
  page: Page, caption: string, columnNames: string[],
): Promise<void> {
  const result = await page.evaluate(({ cap, names }) => {
    const w = window as unknown as { DG?: any; grok?: any };
    if (!w.DG?.Dialog?.getOpenDialogs) return { ok: false, reason: 'DG.Dialog.getOpenDialogs missing' };
    const dialogs: any[] = w.DG.Dialog.getOpenDialogs();
    if (!dialogs.length) return { ok: false, reason: 'no open dialogs' };
    const dlg = dialogs[dialogs.length - 1];
    // `Dialog.inputs` returns raw Dart handles in the current wrapper — re-wrap via DG.toJs so
    // `.caption` and `.value` go through the JS InputBase API.
    const rawInputs: any[] = dlg.inputs ?? [];
    const inputs: any[] = rawInputs.map((d: any) => (w.DG.toJs ? w.DG.toJs(d) : d));
    const captions = inputs.map((i) => String(i?.caption ?? ''));
    const want = cap.toLowerCase();
    let input = inputs.find((i) => String(i?.caption ?? '').toLowerCase() === want);
    if (!input) input = inputs.find((i) => String(i?.caption ?? '').toLowerCase().includes(want));
    if (!input) return { ok: false, reason: `caption "${cap}" not in [${captions.join(', ')}]` };
    const df = w.grok.shell.tv.dataFrame;
    input.value = names.map((n: string) => df.col(n)).filter((c: any) => c != null);
    return { ok: true };
  }, { cap: caption, names: columnNames });
  if (!result.ok) throw new Error(`setDialogColumnListInput: ${result.reason}`);
  await page.waitForTimeout(500);
}

/**
 * Read all numeric column names from the current TableView's DataFrame.
 * Useful when building "all features except X" lists for analysis dialogs.
 */
export async function currentNumericColumnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    const df = g?.shell?.tv?.dataFrame;
    if (!df) return [];
    const out: string[] = [];
    for (let i = 0; i < df.columns.length; i++) {
      const c = df.columns.byIndex(i);
      if (c.type !== 'string') out.push(c.name as string);
    }
    return out;
  });
}

/** Replace the value of a numeric/string input host, confirming with Tab. */
export async function setInputValue(page: Page, hostName: string, value: string): Promise<void> {
  const ed = page.locator(`${inputHost(hostName)} .ui-input-editor`).first();
  await ed.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(value);
  await page.keyboard.press('Tab');
  await page.waitForTimeout(300);
}

/** Toggle a boolean input host on (idempotent — leaves it on if it was already on). */
export async function setBoolInputOn(page: Page, hostName: string): Promise<void> {
  const host = page.locator(inputHost(hostName)).first();
  await host.waitFor({ timeout: 10_000 });
  const isOn = await host.locator('input[type="checkbox"]').first().isChecked().catch(() => false);
  if (!isOn) await host.locator('.ui-input-editor, input[type="checkbox"]').first().click();
  await page.waitForTimeout(200);
}

/**
 * Wait in-browser until every name in `names` is present among the current TableView's
 * columns. Uses `page.waitForFunction` (polls inside the page, rejects only once on
 * timeout) so intermediate "not ready yet" attempts are never logged as report errors —
 * unlike `expect.poll`, which records each failed retry.
 */
export async function waitForColumns(page: Page, names: string[], timeoutMs = 120_000): Promise<void> {
  await page.waitForFunction((wanted: string[]) => {
    const df = (window as any).grok?.shell?.tv?.dataFrame;
    if (!df) return false;
    const have = new Set<string>(df.columns.names() as string[]);
    return wanted.every((n) => have.has(n));
  }, names, { timeout: timeoutMs });
}

/** Return current TableView dataframe column names (small JS-eval — read-only diagnostics). */
export async function currentColumnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    const df = g?.shell?.tv?.dataFrame;
    return df ? (df.columns.names() as string[]) : [];
  });
}

/** Names of currently-attached viewers on the current TableView. */
export async function currentViewerTypes(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    const tv = g?.shell?.tv;
    if (!tv) return [];
    return Array.from(tv.viewers).map((v: any) => v?.type ?? '');
  });
}

/** Visible tab labels rendered inside the page's `d4-tab-host` containers. */
export async function visibleTabLabels(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
      .map((el) => el.textContent?.trim() ?? '')
      .filter((s) => s.length > 0),
  );
}

/** True if the currently-open dialog's primary action button is enabled. */
export async function isPrimaryEnabled(page: Page, candidates: string[] = ['OK', 'Run', 'RUN']): Promise<boolean> {
  for (const c of candidates) {
    const btn = topDialog(page).locator(`[name="button-${c}"]`).first();
    if (await btn.count() > 0) return btn.isEnabled().catch(() => false);
  }
  return false;
}

/** Reset shell between scenarios — close popups, navigate home. */
export async function resetShell(page: Page): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    g?.shell?.closeAll?.();
  }).catch(() => {});
  // No navigation: closeAll() resets the shared page's views/tables in-place, so
  // the next test's openDemoCsv() starts clean without a full reload.
  await page.waitForTimeout(200);
}

/**
 * Open the Train Model view via Top Menu > ML > Models > Train Model... and wait until the
 * platform reports a `PredictiveModel` view in `grok.shell.views`.
 */
export async function openTrainModelView(page: Page): Promise<void> {
  await clickTopMenuLeaf(page, 'div-ML---Models---Train-Model...');
  await page.waitForFunction(() => Array.from((window as any).grok.shell.views)
    .some((v: any) => v?.type === 'PredictiveModel'), undefined, { timeout: 30_000 });
  // Predict input host is the first form field to mount — its presence confirms the form
  // is interactable, not just the view-handle.
  await page.locator(inputHost('Predict')).waitFor({ timeout: 15_000 });
}

/**
 * Type a column name into the Train Model "Predict" picker and confirm with Enter.
 * The picker is a custom column-search popup mounted on `mousedown` of the input editor —
 * `click()` alone races against its focus handlers, so a synthetic `mousedown` is used.
 */
export async function setPredictColumn(page: Page, columnName: string): Promise<void> {
  await page.evaluate(() => {
    const ed = document.querySelector('[name="input-host-Predict"] .ui-input-editor') as HTMLElement | null;
    ed?.dispatchEvent(new MouseEvent('mousedown', { bubbles: true, button: 0 }));
  });
  await page.waitForTimeout(400);
  await page.keyboard.type(columnName);
  await page.keyboard.press('Enter');
  await page.waitForTimeout(400);
}

/**
 * Train an EDA model directly via the function registry, returning a small result envelope
 * useful for assertions. Used as the canonical UI-only-impossible fallback for MLMethods
 * scenarios — see run notes under `public/packages/UsageAnalysis/files/TestTrack/EDA/MLMethods/*-run.md`
 * for the canvas-based features picker and missing Train button issues.
 *
 * Walks `grok.shell.views` to find the first `TableView` (rather than relying on
 * `shell.tv`, which points to the currently-focused view — and that becomes the
 * `PredictiveModel` view once Train Model is opened, dropping the dataset reference).
 */
export async function trainEdaModelViaApi(
  page: Page,
  funcName: string,
  predict: string,
  opts: { numericOnly?: boolean; extraParams?: Record<string, unknown> } = {},
): Promise<{ ok: boolean; error?: string }> {
  return page.evaluate(async ({ fn, target, params, numericOnly }) => {
    try {
      const g = (window as unknown as { grok: any }).grok;
      const tableView = Array.from(g.shell.views).find((v: any) => v?.type === 'TableView') as any;
      if (!tableView?.dataFrame) return { ok: false, error: 'no TableView with a dataFrame in shell.views' };
      const sourceDf = tableView.dataFrame;
      // The predict column always comes from the source frame: with numericOnly a string
      // target (e.g. Species) is dropped from the clone, so reading it off the clone would
      // yield null. Softmax requires numeric-only features but a string target.
      const predictColumn = sourceDf.col(target);
      if (!predictColumn) return { ok: false, error: `predict column "${target}" not found` };
      let df = sourceDf;
      if (numericOnly) {
        const names: string[] = [];
        for (let i = 0; i < df.columns.length; i++) {
          const c = df.columns.byIndex(i);
          if (c.type !== 'string' && c.name !== target) names.push(c.name);
        }
        df = df.clone(null, names);
      }
      const result = await g.functions.call(fn, { df, predictColumn, ...(params ?? {}) });
      return { ok: result != null };
    } catch (e: any) {
      return { ok: false, error: e?.message ?? String(e) };
    }
  }, { fn: funcName, target: predict, params: opts.extraParams ?? null, numericOnly: !!opts.numericOnly });
}

/** Read visible error balloons (empty array = none). */
export async function visibleErrorBalloons(page: Page): Promise<string[]> {
  return page.evaluate(() => Array.from(document.querySelectorAll(
    '.d4-balloon-error, .grok-balloon-error, [class*="balloon"][class*="error"]'))
    .filter((b) => (b as HTMLElement).getBoundingClientRect().width > 0)
    .map((b) => (b as HTMLElement).textContent?.trim() ?? '')
    .filter((s) => s.length > 0));
}

function escapeRegex(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
