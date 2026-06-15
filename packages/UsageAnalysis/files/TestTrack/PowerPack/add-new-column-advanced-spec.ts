/* ---
sub_features_covered: [powerpack.dialogs, powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func, powerpack.dialogs.prepare-add-column-call, powerpack.formula.is-formula-column]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openTableFromFile, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';
test.use(specTestOptions);
test('PowerPack: Add New Column — multi-source datasync persistence + formula recalc on rename (GROK-17109)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;
  const stamp = Date.now();
  const projectName = `AutoTest-AddNewColAdvanced-${stamp}`;
  let projectId: string | null = null;
  let tableInfoId: string | null = null;
  await loginToDatagrok(page);
  await page.evaluate(() => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
  });
  await page.waitForTimeout(500);
  try {
    // Setup: open demog.csv with datasync provenance so the GROK-17109 invariant can be tested.
    await softStep('Setup: open System:DemoFiles/demog.csv with datasync provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);
      // Without wired provenance, save-with-datasync degrades to snapshot-only and GROK-17109 can't be tested.
      await assertProvenanceScript(page, 'files', opened.script);
      const cols = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        return df ? df.columns.names() : [];
      });
      expect(cols).toContain('WEIGHT');
    });
    await softStep('Step 1: open Add New Column dialog via toolbar icon (first time)', async () => {
      const icon = page.locator('[name="icon-add-new-column"]').first();
      await icon.waitFor({timeout: 30_000, state: 'visible'});
      await icon.click({timeout: 10_000});
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await dlg.waitFor({timeout: 30_000});
      await expect(dlg).toBeVisible();
    });
    await softStep('Step 2: add Weight2 = ${WEIGHT}+100 via dialog UI; verify column added', async () => {
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await page.evaluate(() => {
        const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
        if (!input) throw new Error('Name input not found');
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(input, 'Weight2');
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
      });
      await page.waitForTimeout(150);
      const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
      await cm.waitFor({timeout: 15_000, state: 'visible'});
      await cm.click();
      await page.waitForTimeout(200);
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(100);
      let composed: {ok: boolean; doc?: string} = {ok: false};
      for (let i = 0; i < 10; i++) {
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          if (!view) return {ok: false};
          view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: '${WEIGHT} + 100'}});
          return {ok: true, doc: view.state.doc.toString()};
        });
        if (composed.ok) break;
        await page.waitForTimeout(200);
      }
      if (!composed.ok) {
        await cm.click();
        await page.keyboard.press('Control+A');
        await page.keyboard.press('Delete');
        await page.waitForTimeout(100);
        await page.keyboard.type('${WEIGHT} + 100', {delay: 30});
        await page.waitForTimeout(200);
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          const doc = view ? view.state.doc.toString() : (cmDiv.innerText || '');
          return {ok: true, doc};
        });
      }
      expect(composed.ok).toBe(true);
      expect(composed.doc).toContain('${WEIGHT}');
      expect(composed.doc).toContain('+ 100');
      await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
      let added = false;
      for (let i = 0; i < 40; i++) {
        added = await page.evaluate(() => {
          const df = (window as any).grok.shell.tv?.dataFrame;
          return df ? df.columns.names().includes('Weight2') : false;
        });
        if (added) break;
        await page.waitForTimeout(250);
      }
      expect(added).toBe(true);
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w2 = df.col('Weight2');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w2v = w2.get(i);
          if (wv !== null && w2v !== null && Number.isFinite(wv) && Number.isFinite(w2v))
            return {wv, w2v, diff: w2v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(100, 1);
    });
    await softStep('Step 3: reopen Add New Column dialog via toolbar icon (second time)', async () => {
      await page.locator('.d4-dialog').first()
        .waitFor({state: 'detached', timeout: 5_000}).catch(() => {});
      const icon = page.locator('[name="icon-add-new-column"]').first();
      await icon.waitFor({timeout: 15_000, state: 'visible'});
      await icon.click({timeout: 10_000});
      const dlg2 = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await dlg2.waitFor({timeout: 30_000});
      await expect(dlg2).toBeVisible();
    });
    await softStep('Step 4: add Weight3 = ${Weight2}+100 referencing Weight2', async () => {
      const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
      await page.evaluate(() => {
        const input = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement | null;
        if (!input) throw new Error('Name input not found');
        const setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value')!.set!;
        setter.call(input, 'Weight3');
        input.dispatchEvent(new Event('input', {bubbles: true}));
        input.dispatchEvent(new Event('change', {bubbles: true}));
      });
      await page.waitForTimeout(150);
      const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
      await cm.waitFor({timeout: 15_000, state: 'visible'});
      await cm.click();
      await page.waitForTimeout(200);
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(100);
      let composed: {ok: boolean; doc?: string} = {ok: false};
      for (let i = 0; i < 10; i++) {
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          if (!view) return {ok: false};
          view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: '${Weight2} + 100'}});
          return {ok: true, doc: view.state.doc.toString()};
        });
        if (composed.ok) break;
        await page.waitForTimeout(200);
      }
      if (!composed.ok) {
        await cm.click();
        await page.keyboard.press('Control+A');
        await page.keyboard.press('Delete');
        await page.waitForTimeout(100);
        await page.keyboard.type('${Weight2} + 100', {delay: 30});
        await page.waitForTimeout(200);
        composed = await page.evaluate(() => {
          const cmDiv = document.querySelector(
            '.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
          if (!cmDiv) return {ok: false};
          const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
          const doc = view ? view.state.doc.toString() : (cmDiv.innerText || '');
          return {ok: true, doc};
        });
      }
      expect(composed.ok).toBe(true);
      expect(composed.doc).toContain('${Weight2}');
      expect(composed.doc).toContain('+ 100');
      await dlg.locator('[name="button-Add-New-Column---OK"]').first().click();
      let added = false;
      for (let i = 0; i < 40; i++) {
        added = await page.evaluate(() => {
          const df = (window as any).grok.shell.tv?.dataFrame;
          return df ? df.columns.names().includes('Weight3') : false;
        });
        if (added) break;
        await page.waitForTimeout(250);
      }
      expect(added).toBe(true);
      const check = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return null;
        const w = df.col('WEIGHT'); const w3 = df.col('Weight3');
        for (let i = 0; i < Math.min(df.rowCount, 50); i++) {
          const wv = w.get(i); const w3v = w3.get(i);
          if (wv !== null && w3v !== null && Number.isFinite(wv) && Number.isFinite(w3v))
            return {wv, w3v, diff: w3v - wv};
        }
        return null;
      });
      expect(check).not.toBeNull();
      expect(check!.diff).toBeCloseTo(200, 1);
    });
    await softStep('Step 5: rename WEIGHT to BaseWeight via column header context action; edit one cell', async () => {
      const renamed = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df'};
        const col = df.col('WEIGHT');
        if (!col) return {ok: false, why: 'no WEIGHT col'};
        col.name = 'BaseWeight';
        return {ok: true, names: df.columns.names()};
      });
      expect(renamed.ok).toBe(true);
      expect(renamed.names).toContain('BaseWeight');
      expect(renamed.names).not.toContain('WEIGHT');
      await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const col = df.col('BaseWeight');
        const prev = col.get(0);
        col.set(0, (prev ?? 100) + 50);
        df.fireValuesChanged?.();
      });
      await page.waitForTimeout(1000);
      const formula = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2');
        const tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w2v = w2.get(0); const baseV = df.col('BaseWeight').get(0);
        return {tag, w2v, baseV, diff: w2v - baseV};
      });
      expect(formula.tag).toContain('BaseWeight');
      expect(formula.diff).toBeCloseTo(100, 1);
    });
    await softStep('Step 6: save project with datasync provenance preserved', async () => {
      const saved = await saveProjectWithProvenance(page, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoId;
      expect(projectId).toBeTruthy();
      const exists = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        return p != null;
      }, projectId);
      expect(exists).toBe(true);
    });
    await softStep('Step 7: close all views before reopen', async () => {
      await page.evaluate(() => {
        try { (window as any).grok.shell.closeAll(); } catch (_) {}
      });
      await page.waitForTimeout(1000);
      const tableCount = await page.evaluate(() => {
        try { return Number((window as any).grok.shell.tables?.length) || 0; }
        catch { return 0; }
      });
      expect(tableCount).toBe(0);
    });
    await softStep('Step 8: reopen project; verify Weight2 + Weight3 present with formula tags (GROK-17109)', async () => {
      if (!projectId) throw new Error('Step 6 did not produce a projectId');
      const reopen = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        await p.open();
        for (let i = 0; i < 40; i++) {
          const tv = grok.shell.tv;
          if (tv?.dataFrame) break;
          await new Promise((r) => setTimeout(r, 500));
        }
        await new Promise((r) => setTimeout(r, 2000));
        const df = grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df after reopen'};
        const names = df.columns.names();
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        const w2Tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w3Tag = w3?.tags?.get?.('formula') ?? w3?.tags?.get?.('.formula') ?? '';
        return {
          ok: true,
          names,
          hasWeight2: names.includes('Weight2'),
          hasWeight3: names.includes('Weight3'),
          w2Formula: w2Tag,
          w3Formula: w3Tag,
        };
      }, projectId);
      expect(reopen.ok).toBe(true);
      expect(reopen.hasWeight2).toBe(true);
      expect(reopen.hasWeight3).toBe(true);
      expect(reopen.w2Formula.length).toBeGreaterThan(0);
      expect(reopen.w3Formula.length).toBeGreaterThan(0);
    });
    await softStep('Step 9: rename source column post-reopen; verify Weight2 formula updates, Weight3 unaffected', async () => {
      const target = 'BaseWeight2';
      const renamed = await page.evaluate((newName) => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        if (!df) return {ok: false, why: 'no df'};
        const names = df.columns.names();
        const sourceName = names.includes('BaseWeight') ? 'BaseWeight'
                         : names.includes('WEIGHT') ? 'WEIGHT'
                         : null;
        if (!sourceName) return {ok: false, why: 'no source column to rename', names};
        df.col(sourceName).name = newName;
        return {ok: true, sourceName, newNames: df.columns.names()};
      }, target);
      expect(renamed.ok).toBe(true);
      expect(renamed.newNames).toContain(target);
      await page.waitForTimeout(500);
      const tags = await page.evaluate((t) => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        const w2Tag = w2?.tags?.get?.('formula') ?? w2?.tags?.get?.('.formula') ?? '';
        const w3Tag = w3?.tags?.get?.('formula') ?? w3?.tags?.get?.('.formula') ?? '';
        return {w2Tag, w3Tag, target: t};
      }, target);
      expect(tags.w2Tag).toContain(target);
      expect(tags.w3Tag).toContain('Weight2');
    });
    await softStep('Step 10: edit source column values post-reopen; Weight2 + Weight3 recompute', async () => {
      const result = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const sourceName = df.columns.names().find((n: string) =>
          n === 'BaseWeight2' || n === 'BaseWeight' || n === 'WEIGHT') || 'BaseWeight2';
        const src = df.col(sourceName);
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        const preSrc = src.get(0); const preW2 = w2.get(0); const preW3 = w3.get(0);
        const newSrc = (preSrc ?? 100) + 25;
        src.set(0, newSrc);
        df.fireValuesChanged?.();
        return {sourceName, preSrc, preW2, preW3, newSrc};
      });
      await page.waitForTimeout(1000);
      const post = await page.evaluate(() => {
        const df = (window as any).grok.shell.tv?.dataFrame;
        const w2 = df.col('Weight2'); const w3 = df.col('Weight3');
        return {postW2: w2.get(0), postW3: w3.get(0)};
      });
      expect(post.postW2).toBeCloseTo(result.newSrc + 100, 1);
      expect(post.postW3).toBeCloseTo(result.newSrc + 200, 1);
    });
  } finally {
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    await page.evaluate(() => {
      try { (window as any).grok.shell.closeAll(); } catch (_) {}
    }).catch(() => {});
  }
  finishSpec();
});
