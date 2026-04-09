import {test, expect, chromium} from '@playwright/test';

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

test('DiffStudio Scripting: Edit, export JS, run, save, Model Hub', async () => {
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

  // Setup
  await page.evaluate(async () => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = false;
  });

  // Step 1: Open DiffStudio, load Bioreactor, enable Edit, click </>
  await softStep('Step 1: Open DiffStudio, Edit toggle, export JS script', async () => {
    await page!.evaluate(async () => {
      await grok.functions.call('DiffStudio:runDiffStudio');
      await new Promise(r => setTimeout(r, 3000));

      // Load Bioreactor
      const combo = document.querySelector('.diff-studio-ribbon-widget') as HTMLElement;
      if (combo) combo.click();
      await new Promise(r => setTimeout(r, 800));
      let items = document.querySelectorAll('[role="menuitem"]');
      const lib = Array.from(items).find(i => i.textContent?.trim().startsWith('Library'));
      if (lib) (lib as HTMLElement).click();
      await new Promise(r => setTimeout(r, 800));
      items = document.querySelectorAll('[role="menuitem"]');
      const bio = Array.from(items).find(i => i.textContent?.trim() === 'Bioreactor');
      if (bio) (bio as HTMLElement).click();
      await new Promise(r => setTimeout(r, 5000));
    });

    // Click Edit toggle
    await page!.evaluate(async () => {
      const spans = document.querySelectorAll('span.diff-studio-ribbon-text');
      const editSpan = Array.from(spans).find(s => s.textContent?.trim() === 'Edit');
      if (editSpan) (editSpan as HTMLElement).click();
      await new Promise(r => setTimeout(r, 2000));
    });

    // Click </> icon
    await page!.evaluate(async () => {
      const spans = document.querySelectorAll('span.d4-ribbon-name');
      const codeSpan = Array.from(spans).find(s => s.textContent?.trim() === '</>');
      if (codeSpan) (codeSpan as HTMLElement).click();
      await new Promise(r => setTimeout(r, 3000));
    });

    // Verify ScriptView appeared
    const views = await page!.evaluate(() =>
      Array.from(grok.shell.views).map(v => ({name: v.name, type: v.type}))
    );
    expect(views.some(v => v.type === 'ScriptView')).toBe(true);
  });

  // Step 2: Run the script and adjust Final slider
  await softStep('Step 2: Run script and adjust Final input', async () => {
    // Switch to ScriptView and find Run button
    await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const sv = views.find(v => v.type === 'ScriptView');
      if (sv) grok.shell.v = sv;
      await new Promise(r => setTimeout(r, 500));

      // Make play button accessible
      const playIcon = document.querySelector('.d4-ribbon-item i.fa-play');
      if (playIcon) {
        const parent = playIcon.closest('.d4-ribbon-item')!;
        parent.setAttribute('role', 'button');
        parent.setAttribute('aria-label', 'Run script');
      }
    });

    // Click Run
    const runPos = await page!.evaluate(() => {
      const playIcon = document.querySelector('.d4-ribbon-item i.fa-play');
      if (!playIcon) return null;
      const rect = playIcon.getBoundingClientRect();
      return {x: rect.left + rect.width / 2, y: rect.top + rect.height / 2};
    });
    if (runPos) await page!.mouse.click(runPos.x, runPos.y);
    await page!.waitForTimeout(8000);

    // Switch to the js-view-base result
    await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const jsView = views.find(v => v.type === 'js-view-base');
      if (jsView) grok.shell.v = jsView;
      await new Promise(r => setTimeout(r, 2000));
    });

    // Modify Final input
    const result = await page!.evaluate(async () => {
      const finalInput = document.querySelector('input[name="input-Final"]') as HTMLInputElement;
      if (!finalInput) return {error: 'not found'};
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(finalInput, '500');
      finalInput.dispatchEvent(new Event('input', {bubbles: true}));
      finalInput.dispatchEvent(new Event('change', {bubbles: true}));
      finalInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
      return {value: finalInput.value};
    });
    expect(result.value).toBe('500');
  });

  // Step 3: Add //tags: model and save
  await softStep('Step 3: Add tags and save script', async () => {
    // Switch to ScriptView
    await page!.evaluate(async () => {
      const views = Array.from(grok.shell.views);
      const sv = views.find(v => v.type === 'ScriptView');
      if (sv) grok.shell.v = sv;
      await new Promise(r => setTimeout(r, 1000));
    });

    // Add //tags: model via CodeMirror
    const tagResult = await page!.evaluate(() => {
      const cm = document.querySelector('.CodeMirror') as any;
      if (!cm?.CodeMirror) return {error: 'CodeMirror not found'};
      const content = cm.CodeMirror.getValue();
      if (!content.includes('//tags: model')) {
        const lines = content.split('\n');
        const nameIdx = lines.findIndex((l: string) => l.startsWith('//name:'));
        const tagsIdx = lines.findIndex((l: string) => l.startsWith('//tags:'));
        if (tagsIdx >= 0) {
          if (!lines[tagsIdx].includes('model'))
            lines[tagsIdx] = lines[tagsIdx].trimEnd() + ', model';
        } else if (nameIdx >= 0) {
          lines.splice(nameIdx + 1, 0, '//tags: model');
        }
        cm.CodeMirror.setValue(lines.join('\n'));
      }
      return {hasTag: cm.CodeMirror.getValue().includes('//tags: model')};
    });
    expect(tagResult.hasTag).toBe(true);

    // Click SAVE
    const savePos = await page!.evaluate(() => {
      const btns = document.querySelectorAll('button');
      const saveBtn = Array.from(btns).find(b => b.textContent?.trim() === 'SAVE');
      if (!saveBtn) return null;
      const rect = saveBtn.getBoundingClientRect();
      return {x: rect.left + rect.width / 2, y: rect.top + rect.height / 2};
    });
    if (savePos) await page!.mouse.click(savePos.x, savePos.y);
    await page!.waitForTimeout(3000);
  });

  // Step 4: Access Model in Model Hub
  await softStep('Step 4: Open Model Hub and find Bioreactor', async () => {
    const viewName = await page!.evaluate(async () => {
      const v = await grok.functions.call('Compute2:modelCatalog');
      grok.shell.v = v;
      await new Promise(r => setTimeout(r, 5000));
      return grok.shell.v?.name;
    });
    expect(viewName).toBe('Model Hub');

    // Verify Bioreactor card exists
    const found = await page!.evaluate(() => {
      const labels = document.querySelectorAll('span.d4-link-label');
      return Array.from(labels).some(l => l.textContent?.trim() === 'Bioreactor');
    });
    expect(found).toBe(true);
  });

  // Step 5: Open model from Hub and interact
  await softStep('Step 5: Open Bioreactor from Hub and adjust Final', async () => {
    await page!.evaluate(async () => {
      const labels = document.querySelectorAll('span.d4-link-label');
      for (const el of labels) {
        if (el.textContent?.trim() === 'Bioreactor') {
          (el as HTMLElement).click();
          break;
        }
      }
      await new Promise(r => setTimeout(r, 1000));
      const obj = grok.shell.o;
      if (obj && typeof obj.apply === 'function') await obj.apply();
      await new Promise(r => setTimeout(r, 5000));

      const views = Array.from(grok.shell.views);
      const jsViews = views.filter(v => v.type === 'js-view-base');
      if (jsViews.length > 0) grok.shell.v = jsViews[jsViews.length - 1];
      await new Promise(r => setTimeout(r, 2000));
    });

    // Modify Final input
    const result = await page!.evaluate(async () => {
      const finalInput = document.querySelector('input[name="input-Final"]') as HTMLInputElement;
      if (!finalInput) return {error: 'not found'};
      const nativeSetter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
      nativeSetter.call(finalInput, '600');
      finalInput.dispatchEvent(new Event('input', {bubbles: true}));
      finalInput.dispatchEvent(new Event('change', {bubbles: true}));
      finalInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', bubbles: true}));
      await new Promise(r => setTimeout(r, 2000));
      const canvases = document.querySelectorAll('canvas').length;
      return {value: finalInput.value, canvases};
    });
    expect(result.value).toBe('600');
    expect(result.canvases).toBeGreaterThanOrEqual(1);
  });

  // Cleanup: delete saved script and close all
  await page.evaluate(async () => {
    try {
      const scripts = await grok.dapi.scripts.filter('name = "Bioreactor"').list();
      const myScripts = scripts.filter((s: any) => s.author?.login === 'claude');
      for (const s of myScripts) await grok.dapi.scripts.delete(s);
    } catch(e) {}
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
