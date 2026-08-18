import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

interface Probe {
  query: string;
  branchHint: string;
  isCanonicalRepro: boolean;
}

const DISPATCH_PROBES: Probe[] = [

  {query: 'QA', branchHint: 'short-query; same shape that triggered GROK-18656', isCanonicalRepro: true},

  {query: 'a', branchHint: 'single-char short-query — most aggressive null-safety surface', isCanonicalRepro: false},
  {query: 'new', branchHint: 'common word — function-call / widget dispatch branches', isCanonicalRepro: false},
  {query: 'user', branchHint: 'username-style — regex-entity (users) + function-call branches', isCanonicalRepro: false},
  {query: 'Project[0-9]+', branchHint: 'regex-pattern entity — regex-entity dispatch branch', isCanonicalRepro: false},
  {query: '1+1', branchHint: 'function-evaluation expression — JS-eval / function-evaluation branch', isCanonicalRepro: false},
];

test('PowerPack: Power Search Enter-key dispatch is null-safe across 8 paths (GROK-18656 regression)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const pageErrors: Array<{message: string; stack?: string}> = [];
  page.on('pageerror', (err) => {
    pageErrors.push({message: err.message, stack: err.stack});
  });
  const reportedErrors: string[] = [];
  page.on('console', (msg) => {
    if (msg.type() !== 'error') return;
    const text = msg.text();
    if (/reported error/i.test(text) || /Uncaught/i.test(text))
      reportedErrors.push(text);
  });

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    try { grok.shell.windows.simpleMode = true; } catch (_) {}
    try { grok.shell.closeAll(); } catch (_) {}
    try {
      const v = await grok.functions.call('PowerPack:_welcomeView', {});
      if (v) grok.shell.addView(v);
    } catch (_) {  }
  });

  const welcomeView = page.locator('.power-pack-welcome-view').first();
  await welcomeView.waitFor({timeout: 60_000, state: 'visible'});

  const searchInput = page.locator('input.power-search-search-everywhere-input').first();
  await searchInput.waitFor({timeout: 30_000, state: 'visible'});

  await page.waitForTimeout(1500); 

  async function clearSearchInput(p: Page): Promise<void> {
    const input = p.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    await p.keyboard.press('Control+A');
    await p.keyboard.press('Delete');
    await p.waitForTimeout(150);
  }

  function errorSnapshot(): {pageErrorsCount: number; reportedErrorsCount: number} {
    return {pageErrorsCount: pageErrors.length, reportedErrorsCount: reportedErrors.length};
  }

  async function readShellLastError(p: Page): Promise<string | null> {
    return p.evaluate(() => {
      const g = (window as any).grok;
      try {
        const v = g?.shell?.lastError;
        return typeof v === 'string' ? v : (v == null ? null : String(v));
      } catch (_) { return null; }
    });
  }

  function assertNoNewErrors(
    label: string,
    baseline: {pageErrorsCount: number; reportedErrorsCount: number},
  ): void {
    const newPageErrors = pageErrors.slice(baseline.pageErrorsCount);
    const newReported = reportedErrors.slice(baseline.reportedErrorsCount);
    if (newPageErrors.length > 0) {
      const summary = newPageErrors.map((e) => `  ${e.message}`).join('\n');
      throw new Error(
        `GROK-18656 invariant violated for probe '${label}': ` +
        `${newPageErrors.length} new uncaught pageerror raised during dispatch:\n${summary}`,
      );
    }
    if (newReported.length > 0) {
      const hasReportedNull = newReported.some((m) => /reported error\s*:\s*null/i.test(m));
      if (hasReportedNull) {
        throw new Error(
          `GROK-18656 invariant violated for probe '${label}': ` +
          `'reported error: null' surfaced on the console:\n  ${newReported.join('\n  ')}`,
        );
      }

    }
  }

  async function dispatchProbe(p: Page, probe: Probe): Promise<void> {
    await clearSearchInput(p);
    const baseline = errorSnapshot();
    const input = p.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    await p.keyboard.type(probe.query, {delay: 40});
    await p.waitForTimeout(900); 
    await p.keyboard.press('Enter');
    await p.waitForTimeout(700); 
    assertNoNewErrors(probe.query, baseline);
    const shellErr = await readShellLastError(p);
    if (shellErr && /reported error\s*:\s*null/i.test(shellErr))
      throw new Error(
        `GROK-18656 invariant violated for probe '${probe.query}': ` +
        `grok.shell.lastError contains 'reported error: null': ${shellErr}`,
      );
  }

  await softStep('Scenario 1: focus + type "QA" + press Enter → dispatch must not throw (GROK-18656 invariant)', async () => {
    const probe = DISPATCH_PROBES[0]; 
    expect(probe.isCanonicalRepro).toBe(true);
    await dispatchProbe(page, probe);
  });

  await softStep('Scenario 1 Step 5: verify a sensible dispatch result is rendered (results, empty-state, or graceful balloon)', async () => {

    const searchHost = page.locator('.power-pack-search-host').first();
    await searchHost.waitFor({timeout: 10_000, state: 'visible'}).catch(() => {});
    const hostVisible = await searchHost.isVisible({timeout: 1000}).catch(() => false);
    const hostChildCount = hostVisible
      ? await page.evaluate(() => {
        const el = document.querySelector('.power-pack-search-host');
        return el ? el.children.length : 0;
      })
      : 0;
    const balloonCount = await page.evaluate(() => {
      return document.querySelectorAll('.grok-balloon, .d4-balloon').length;
    });
    const sensible = (hostVisible && hostChildCount > 0) ||
                     (hostVisible && hostChildCount === 0) ||
                     (balloonCount > 0);
    expect(sensible, `expected a sensible dispatch result for 'QA' — search-host visible: ${hostVisible}, ` +
      `host child count: ${hostChildCount}, balloon count: ${balloonCount}`).toBe(true);
  });

  for (let i = 1; i < DISPATCH_PROBES.length; i++) {
    const probe = DISPATCH_PROBES[i];
    await softStep(
      `Scenario 2 probe ${i + 1}/${DISPATCH_PROBES.length}: ` +
      `focus + type ${JSON.stringify(probe.query)} + Enter → dispatch must not throw ` +
      `(${probe.branchHint})`,
      async () => {
        await dispatchProbe(page, probe);
      },
    );
  }

  await softStep('Scenario 3 Step 1-2: focus search input + type "dem" to surface suggestion menu', async () => {
    await clearSearchInput(page);
    const input = page.locator('input.power-search-search-everywhere-input').first();
    await input.click({timeout: 10_000});
    await page.keyboard.type('dem', {delay: 50});
    await page.waitForTimeout(900);
  });

  await softStep('Scenario 3 Step 3: press ArrowDown to highlight first suggestion (if any present)', async () => {
    const baseline = errorSnapshot();
    const itemsBeforeCount = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return 0;
      return root.querySelectorAll('.d4-menu-item').length;
    });
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(300);
    assertNoNewErrors('ArrowDown-after-dem', baseline);
    if (itemsBeforeCount > 0) {
      const highlightedCount = await page.evaluate(() => {
        const root = document.querySelector('.power-pack-welcome-view');
        if (!root) return 0;
        return root.querySelectorAll('.d4-menu-item-hover').length;
      });
      expect(highlightedCount).toBeLessThanOrEqual(itemsBeforeCount);
    }
  });

  await softStep('Scenario 3 Step 4: press ArrowDown again (or ArrowUp) — navigation in both directions does not throw', async () => {
    const baseline = errorSnapshot();
    await page.keyboard.press('ArrowDown');
    await page.waitForTimeout(200);
    await page.keyboard.press('ArrowUp');
    await page.waitForTimeout(200);
    assertNoNewErrors('ArrowDown+ArrowUp-after-dem', baseline);
  });

  await softStep('Scenario 3 Step 5: press Enter on highlighted suggestion (if any) — activation does not throw, GROK-18656 fix preserves this path', async () => {
    const baseline = errorSnapshot();
    const hadHighlighted = await page.evaluate(() => {
      const root = document.querySelector('.power-pack-welcome-view');
      if (!root) return false;
      return root.querySelector('.d4-menu-item-hover') != null;
    });
    await page.keyboard.press('Enter');
    await page.waitForTimeout(800);
    assertNoNewErrors('Enter-on-highlighted-suggestion', baseline);
    void hadHighlighted; 
  });

  try {
    await clearSearchInput(page);
  } catch (_) {  }
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) {  }
  }).catch(() => {});

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    const errCount = pageErrors.length;
    const reportedCount = reportedErrors.length;
    throw new Error(
      `${stepErrors.length} step(s) failed (pageerror count: ${errCount}; ` +
      `reported-error count: ${reportedCount}):\n${summary}`,
    );
  }
});
