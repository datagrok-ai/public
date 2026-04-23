import {test, expect} from '@playwright/test';

test.use({
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 60_000,
});

const baseUrl = process.env.DATAGROK_URL ?? 'http://localhost:8888';
const login = process.env.DATAGROK_LOGIN ?? 'admin';
const password = process.env.DATAGROK_PASSWORD ?? 'admin';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try { await test.step(name, fn); }
  catch (e: any) { stepErrors.push({step: name, error: e?.message ?? String(e)}); }
}

test('Database meta — sticky meta on DB connections/tables/columns', async ({page}) => {
  test.setTimeout(300_000);

  await page.goto(baseUrl);
  const loginInput = page.getByPlaceholder('Login or Email').and(page.locator(':visible'));
  if (await loginInput.isVisible({timeout: 15000}).catch(() => false)) {
    await loginInput.click();
    await page.keyboard.type(login);
    await page.getByPlaceholder('Password').and(page.locator(':visible')).click();
    await page.keyboard.type(password);
    await page.keyboard.press('Enter');
  }
  await page.locator('[name="Browse"]').waitFor({timeout: 120000});

  // Setup phase
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    g.shell.windows.showBrowse = true;
    await new Promise((r) => setTimeout(r, 500));
  });

  // ---- Database metadata display section ----
  await softStep('DB: open Browse > Databases > Postgres > CHEMBL', async () => {
    await page.evaluate(async () => {
      const click = (sel: string) => (document.querySelector(sel) as HTMLElement | null)?.click();
      click('[name="tree-expander-Databases"]');
      await new Promise((r) => setTimeout(r, 400));
      click('[name="tree-expander-Databases---Postgres"]');
      await new Promise((r) => setTimeout(r, 800));
      click('[name="tree-Databases---Postgres---CHEMBL"]');
      await new Promise((r) => setTimeout(r, 2500));
    });
    const selectedName = await page.evaluate(() => (window as any).grok.shell.o?.name);
    expect(selectedName).toBe('Chembl');
  });

  await softStep('DB: Context Panel shows "Database meta" section', async () => {
    const panels = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.grok-entity-prop-panel .d4-accordion-pane-header'))
        .map((h) => (h as HTMLElement).textContent?.trim() ?? '')
    );
    expect(panels.some((p) => /Database\s*meta/i.test(p))).toBe(true);
  });

  await softStep('DB: fill Comment=test, LLM Comment=test, Save', async () => {
    const result = await page.evaluate(async () => {
      const hdr = Array.from(document.querySelectorAll('.grok-entity-prop-panel .d4-accordion-pane-header'))
        .find((h) => /Database\s*meta/i.test(h.textContent || ''));
      if (!hdr) return {ok: false, error: 'no Database meta section'};
      (hdr as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 300));
      const commentInput = document.querySelector('[name="input-Comment"] input, [name="input-host-Comment"] input') as HTMLInputElement | null;
      const llmInput = document.querySelector('[name="input-LLM-Comment"] input, [name="input-host-LLM-Comment"] input, [name="input-llmComment"] input') as HTMLInputElement | null;
      if (!commentInput || !llmInput) return {ok: false, error: 'inputs not found'};
      commentInput.value = 'test'; commentInput.dispatchEvent(new Event('input', {bubbles: true}));
      llmInput.value = 'test'; llmInput.dispatchEvent(new Event('input', {bubbles: true}));
      const save = document.querySelector('.grok-entity-prop-panel [name="button-Save"]') as HTMLElement | null;
      save?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return {ok: true};
    });
    expect((result as any).ok).toBe(true);
  });

  await softStep('DB: reload and verify saved values persist', async () => {
    await page.reload();
    await page.locator('[name="Browse"]').waitFor({timeout: 60000});
    await page.evaluate(async () => {
      (window as any).grok.shell.windows.showBrowse = true;
      const click = (sel: string) => (document.querySelector(sel) as HTMLElement | null)?.click();
      click('[name="tree-expander-Databases"]');
      await new Promise((r) => setTimeout(r, 400));
      click('[name="tree-expander-Databases---Postgres"]');
      await new Promise((r) => setTimeout(r, 800));
      click('[name="tree-Databases---Postgres---CHEMBL"]');
      await new Promise((r) => setTimeout(r, 2500));
    });
    const values = await page.evaluate(() => {
      const c = document.querySelector('[name="input-Comment"] input, [name="input-host-Comment"] input') as HTMLInputElement | null;
      const l = document.querySelector('[name="input-LLM-Comment"] input, [name="input-host-LLM-Comment"] input, [name="input-llmComment"] input') as HTMLInputElement | null;
      return {comment: c?.value, llm: l?.value};
    });
    expect(values.comment).toBe('test');
    expect(values.llm).toBe('test');
  });

  await softStep('DB: clear values, save, reload, verify empty', async () => {
    await page.evaluate(async () => {
      const c = document.querySelector('[name="input-Comment"] input, [name="input-host-Comment"] input') as HTMLInputElement | null;
      const l = document.querySelector('[name="input-LLM-Comment"] input, [name="input-host-LLM-Comment"] input, [name="input-llmComment"] input') as HTMLInputElement | null;
      if (c) { c.value = ''; c.dispatchEvent(new Event('input', {bubbles: true})); }
      if (l) { l.value = ''; l.dispatchEvent(new Event('input', {bubbles: true})); }
      const save = document.querySelector('.grok-entity-prop-panel [name="button-Save"]') as HTMLElement | null;
      save?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.reload();
    await page.locator('[name="Browse"]').waitFor({timeout: 60000});
    const values = await page.evaluate(() => {
      const c = document.querySelector('[name="input-Comment"] input, [name="input-host-Comment"] input') as HTMLInputElement | null;
      const l = document.querySelector('[name="input-LLM-Comment"] input, [name="input-host-LLM-Comment"] input, [name="input-llmComment"] input') as HTMLInputElement | null;
      return {comment: c?.value ?? '', llm: l?.value ?? ''};
    });
    expect(values.comment).toBe('');
    expect(values.llm).toBe('');
  });

  // ---- Table metadata display section ----
  await softStep('Table: open NorthwindTest > Schemas > Public', async () => {
    const hasSchemasNode = await page.evaluate(async () => {
      const nw = document.querySelector('[name="tree-Databases---Postgres---NorthwindTest"]') as HTMLElement | null;
      nw?.click();
      await new Promise((r) => setTimeout(r, 2500));
      return !!document.querySelector('[name*="NorthwindTest---Schemas"]');
    });
    expect(hasSchemasNode).toBe(true);
  });

  await softStep('Table: Context Panel shows "Database meta" section', async () => {
    const hasSection = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.grok-entity-prop-panel .d4-accordion-pane-header'))
        .some((h) => /Database\s*meta/i.test(h.textContent || ''))
    );
    expect(hasSection).toBe(true);
  });

  // ---- Column metadata display section ----
  await softStep('Column: open categories > categoryid', async () => {
    await page.goto(`${baseUrl}/dbtable/Postgres.PostgresTest.public.categories?browse=db`);
    await page.waitForTimeout(4000);
    const viewName = await page.evaluate(() => (window as any).grok.shell.v?.name);
    expect(viewName).not.toBe('Home');
  });

  await softStep('Column: Context Panel shows "Database meta" section', async () => {
    const hasSection = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.grok-entity-prop-panel .d4-accordion-pane-header'))
        .some((h) => /Database\s*meta/i.test(h.textContent || ''))
    );
    expect(hasSection).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Scenario failures:\n' + stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
});
