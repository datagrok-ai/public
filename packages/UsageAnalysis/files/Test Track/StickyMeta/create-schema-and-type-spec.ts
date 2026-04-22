import {test, expect, chromium} from '@playwright/test';

test.use({
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('StickyMeta Create Schema: Navigate to Schemas, verify TestSchema1', async () => {
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  await page.evaluate(() => {
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
  });

  // Step 1: Navigate to Sticky Meta Schemas
  await softStep('Step 1: Navigate to Browse > Platform > Sticky Meta > Schemas', async () => {
    await page!.goto(`${baseUrl}/meta/schemas`, {waitUntil: 'networkidle', timeout: 15000});
    await page!.waitForTimeout(3000);

    const pageTitle = await page!.evaluate(() => document.title);
    expect(pageTitle).toContain('Schemas');
  });

  // Step 2: Verify TestSchema1 exists
  await softStep('Step 2: Verify TestSchema1 in schemas list', async () => {
    const hasTestSchema = await page!.evaluate(() => {
      const allEls = document.querySelectorAll('*');
      return Array.from(allEls).some(el =>
        el.textContent?.trim() === 'TestSchema1' && el.children.length === 0
      );
    });
    expect(hasTestSchema).toBe(true);
  });

  // Cleanup
  await page.evaluate(() => grok.shell.closeAll());

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
