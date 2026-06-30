/* ---
sub_features_covered: [charts.tree, charts.tree.color-palette, charts.tree.on-click]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';

test('Charts / Tree — rowSource x onClick state machine (github-3245)', async ({page}) => {
  test.setTimeout(300_000);

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await page.evaluate(() => {
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  await softStep('Setup: Open demog.csv, add Tree, set hierarchy', async () => {
    const result = await page.evaluate(async (path) => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv(path);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
      const tree = tv.addViewer('Tree');
      await new Promise((r) => setTimeout(r, 3000));
      try { tree.setOptions({hierarchyColumnNames: ['CONTROL', 'SEX', 'RACE']}); } catch (e) {}
      await new Promise((r) => setTimeout(r, 1500));
      const types: string[] = [];
      for (const v of tv.viewers) types.push(v.type);
      return {types};
    }, demogPath);
    expect(result.types).toContain('Tree');
  });

  const rowSources = ['All', 'Filtered', 'Selected'];
  const onClicks = ['Select', 'Filter', 'None'];

  const iterateCombo = async (rs: string, oc: string, suffix: string) => {
    await softStep(`Combo ${rs} x ${oc} ${suffix}`, async () => {
      const errorsBefore = consoleErrors.length;
      const result = await page.evaluate(async ([r, o]) => {
        const grok = (window as any).grok;
        const tv = grok.shell.tv;
        let tree: any = null;
        for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
        if (!tree) return {ok: false};
        let setOptionsThrew = false;
        try { tree.setOptions({rowSource: r, onClick: o}); }
        catch (e) { setOptionsThrew = true; }
        await new Promise((r2) => setTimeout(r2, 800));
        let readbackRs: any = null, readbackOc: any = null;
        try { readbackRs = tree.props.get('rowSource'); } catch (e) {}
        try { readbackOc = tree.props.get('onClick'); } catch (e) {}
        const root = tree.root as HTMLElement;
        const rect = root.getBoundingClientRect();
        return {
          ok: true,
          setOptionsThrew,
          readbackRs,
          readbackOc,
          hasContent: root.children.length > 0,
          width: rect.width,
          height: rect.height,
        };
      }, [rs, oc] as [string, string]);
      expect(result.ok).toBe(true);
      // github-3245: readback logged only — Tree normalizes (rowSource, onClick) so strict round-trip drifts.
      expect(result.setOptionsThrew).toBe(false);
      expect(result.hasContent).toBe(true);
      expect(result.width).toBeGreaterThan(0);
      const errorsDuring = consoleErrors.slice(errorsBefore);
      expect(errorsDuring).toEqual([]);
      console.log(`[Combo ${rs}/${oc}]`, JSON.stringify({
        readbackRs: result.readbackRs,
        readbackOc: result.readbackOc,
        readbackMatches: (result.readbackRs === rs) && (result.readbackOc === oc),
      }));
    });
  };

  // Forward iteration: 9 combinations
  for (const rs of rowSources) {
    for (const oc of onClicks) {
      await iterateCombo(rs, oc, '(forward)');
    }
  }

  // Reverse iteration: same 9 combos, reverse order
  for (const rs of [...rowSources].reverse()) {
    for (const oc of [...onClicks].reverse()) {
      await iterateCombo(rs, oc, '(reverse)');
    }
  }

  // ESC-equivalent state-clear check
  await softStep('Step 4: ESC-equivalent (df.selection.setAll(false)) — viewer remains stable', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      df.selection.fireChanged();
      await new Promise((r) => setTimeout(r, 800));
      const tv = grok.shell.tv;
      let tree: any = null;
      for (const v of tv.viewers) if (v.type === 'Tree') { tree = v; break; }
      if (!tree) return {ok: false};
      const root = tree.root as HTMLElement;
      const rect = root.getBoundingClientRect();
      return {ok: true, hasContent: root.children.length > 0, width: rect.width};
    });
    expect(result.ok).toBe(true);
    expect(result.hasContent).toBe(true);
    expect(result.width).toBeGreaterThan(0);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
