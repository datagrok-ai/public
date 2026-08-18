import { test as base, expect, type Page } from '@playwright/test';

export const BASE = process.env.DATAGROK_URL!;

export const test = base.extend<{}, { sharedPage: Page }>({
  sharedPage: [async ({ browser }, use) => {
    const context = await browser.newContext({ storageState: 'e2e/.auth.json' });
    const page = await context.newPage();
    await page.goto(BASE, { waitUntil: 'domcontentloaded', timeout: 120_000 });
    await page.waitForFunction(() => document.querySelector('.grok-preloader') == null,
      undefined, { timeout: 120_000 });
    await page.locator('.d4-ribbon').first().waitFor({ timeout: 60_000 });
    await use(page);
    await context.close();
  }, { scope: 'worker' }],

  page: async ({ sharedPage }, use) => {
    await use(sharedPage);
  },
});

export { expect };

export const inputHost = (safeName: string): string => `[name="input-host-${safeName}"]`;

export const inputEditor = (safeName: string): string => `${inputHost(safeName)} .ui-input-editor`;

export async function openDemoCsv(page: Page, fileName: string, uiTimeoutMs = 25_000): Promise<void> {
  await page.evaluate(async (name) => {
    const g = (window as unknown as { grok: any }).grok;
    const df = await g.dapi.files.readCsv(`System:DemoFiles/${name}`);
    g.shell.addTableView(df);
  }, fileName);
  await waitForCurrentTableView(page, uiTimeoutMs);
}

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

export async function waitForDialog(page: Page, title: string, timeoutMs = 15_000): Promise<void> {
  const re = new RegExp(escapeRegex(title), 'i');
  await page.locator('.d4-dialog .d4-dialog-title', { hasText: re }).first().waitFor({ timeout: timeoutMs });
}

export function topDialog(page: Page) {
  return page.locator('.d4-dialog').last();
}

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

export async function selectAllColumnsInPicker(page: Page, hostName: string): Promise<void> {
  await page.locator(`${inputHost(hostName)} .ui-input-editor`).first().click();

  await page.locator('.d4-dialog [name="label-All"]').last().waitFor({ timeout: 10_000 });
  await page.locator('.d4-dialog [name="label-All"]').last().click();
  await page.locator('.d4-dialog [name="button-OK"]').last().click();
  await page.waitForTimeout(500);
}

export async function setDialogColumnListInput(
  page: Page, caption: string, columnNames: string[],
): Promise<void> {
  const result = await page.evaluate(({ cap, names }) => {
    const w = window as unknown as { DG?: any; grok?: any };
    if (!w.DG?.Dialog?.getOpenDialogs) return { ok: false, reason: 'DG.Dialog.getOpenDialogs missing' };
    const dialogs: any[] = w.DG.Dialog.getOpenDialogs();
    if (!dialogs.length) return { ok: false, reason: 'no open dialogs' };
    const dlg = dialogs[dialogs.length - 1];

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

export async function setInputValue(page: Page, hostName: string, value: string): Promise<void> {
  const ed = page.locator(`${inputHost(hostName)} .ui-input-editor`).first();
  await ed.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(value);
  await page.keyboard.press('Tab');
  await page.waitForTimeout(300);
}

export async function setBoolInputOn(page: Page, hostName: string): Promise<void> {
  const host = page.locator(inputHost(hostName)).first();
  await host.waitFor({ timeout: 10_000 });
  const isOn = await host.locator('input[type="checkbox"]').first().isChecked().catch(() => false);
  if (!isOn) await host.locator('.ui-input-editor, input[type="checkbox"]').first().click();
  await page.waitForTimeout(200);
}

export async function waitForColumns(page: Page, names: string[], timeoutMs = 120_000): Promise<void> {
  await page.waitForFunction((wanted: string[]) => {
    const df = (window as any).grok?.shell?.tv?.dataFrame;
    if (!df) return false;
    const have = new Set<string>(df.columns.names() as string[]);
    return wanted.every((n) => have.has(n));
  }, names, { timeout: timeoutMs });
}

export async function currentColumnNames(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    const df = g?.shell?.tv?.dataFrame;
    return df ? (df.columns.names() as string[]) : [];
  });
}

export async function currentViewerTypes(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    const tv = g?.shell?.tv;
    if (!tv) return [];
    return Array.from(tv.viewers).map((v: any) => v?.type ?? '');
  });
}

export async function visibleTabLabels(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
      .map((el) => el.textContent?.trim() ?? '')
      .filter((s) => s.length > 0),
  );
}

export async function isPrimaryEnabled(page: Page, candidates: string[] = ['OK', 'Run', 'RUN']): Promise<boolean> {
  for (const c of candidates) {
    const btn = topDialog(page).locator(`[name="button-${c}"]`).first();
    if (await btn.count() > 0) return btn.isEnabled().catch(() => false);
  }
  return false;
}

export async function resetShell(page: Page): Promise<void> {
  await page.keyboard.press('Escape').catch(() => {});
  await page.evaluate(() => {
    const g = (window as unknown as { grok?: any }).grok;
    g?.shell?.closeAll?.();
  }).catch(() => {});

  await page.waitForTimeout(200);
}

export async function openTrainModelView(page: Page): Promise<void> {
  await clickTopMenuLeaf(page, 'div-ML---Models---Train-Model...');
  await page.waitForFunction(() => Array.from((window as any).grok.shell.views)
    .some((v: any) => v?.type === 'PredictiveModel'), undefined, { timeout: 30_000 });

  await page.locator(inputHost('Predict')).waitFor({ timeout: 15_000 });
}

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
