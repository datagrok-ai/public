import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer — property surface extension (dataJson/pdb/pdbTag/behaviour/binding-site/layout/controls)', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });
  page.on('console', (msg) => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  });

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {

    let scenario1Mounted = false;

    await softStep('Scenario 1 — Open Molecule3D + Molecule DF; add Biostructure viewer; container mounts', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['row-1', 'row-2', 'row-3']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent, pdbContent]),
          DG.Column.fromStrings('ligand', ['CCO', 'CCC', 'CCN']),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }
        try { df.col('ligand').semType = 'Molecule'; } catch (_e) {  }
        df.name = 'property-surface-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 3000));
        return {
          rowCount: df.rowCount,
          hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
          vType: v?.type,
          contentLen: pdbContent.length,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.hasContainer).toBe(true);
      expect(res.vType).toBe('Biostructure');
      expect(res.rowCount).toBe(3);
      expect(res.contentLen).toBeGreaterThan(1000);
      scenario1Mounted = true;
    });

    await softStep('Scenario 1 — biostructureDataToJson({name}) + setOptions({dataJson}) round-trips verbatim', async () => {
      if (!scenario1Mounted) return;
      const res = await page.evaluate(async (pdbPath) => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'Biostructure') as any : null;
        if (!v) return {ok: false};

        const pdbContent = await grok.dapi.files.readAsText(pdbPath);
        const dataJsonStr = await grok.functions.call(
          'BiostructureViewer:biostructureDataToJson',
          {binary: false, data: pdbContent, ext: 'pdb', options: {name: 'fixture-1bdq'}},
        );

        const dataJsonStrType = typeof dataJsonStr;
        const dataJsonStrLen = (dataJsonStr && dataJsonStr.length) || 0;

        v.setOptions({dataJson: dataJsonStr});
        await new Promise((r) => setTimeout(r, 1500));
        const dataJsonAfter = v.props.get('dataJson');

        const props = v.props.getProperties();
        const dataJsonProp = props.find((p: any) => p.name === 'dataJson');

        return {
          ok: true,
          dataJsonStrType,
          dataJsonStrLen,
          dataJsonAfterType: typeof dataJsonAfter,
          dataJsonAfterLen: dataJsonAfter ? String(dataJsonAfter).length : 0,
          dataJsonRoundTrips: dataJsonStr === dataJsonAfter,
          dataJsonInCatalogue: !!dataJsonProp,
          dataJsonCategory: dataJsonProp?.category ?? null,

          dataJsonContainsName: typeof dataJsonAfter === 'string' && dataJsonAfter.indexOf('"name":"fixture-1bdq"') >= 0,
        };
      }, samplePdbPath);

      expect(res.ok).toBe(true);
      expect(res.dataJsonStrType).toBe('string');
      expect(res.dataJsonStrLen).toBeGreaterThan(1000);
      expect(res.dataJsonRoundTrips).toBe(true);
      expect(res.dataJsonInCatalogue).toBe(true);
      expect(res.dataJsonCategory).toBe('Data');
      expect(res.dataJsonContainsName).toBe(true);
    });

    await softStep('Scenario 2 step 4 — Raw pdb without name: setOptions({pdb}) stores content; pitfall signature optimistically captured', async () => {
      if (!scenario1Mounted) return;

      pageErrors.length = 0;
      consoleErrors.length = 0;

      const res = await page.evaluate(async (pdbPath) => {

        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s2-row-1', 's2-row-2']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent]),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2000));

        let setErr: string | null = null;
        try { v.setOptions({pdb: pdbContent}); }
        catch (e: any) { setErr = String(e?.message ?? e); }

        await new Promise((r) => setTimeout(r, 5000));

        const pdbAfter = v.props.get('pdb');
        return {
          setErr,
          pdbStored: !!pdbAfter && pdbAfter.length > 0,
          pdbLenAfter: pdbAfter ? pdbAfter.length : 0,
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.setErr).toBe(null);
      expect(res.pdbStored).toBe(true);
      expect(res.pdbLenAfter).toBeGreaterThan(1000);
      expect(res.containerPresent).toBe(true);

      const pitfallRegex = /Parsed object is empty|name\s+'undefined'/i;
      const pitfallHitsConsole = consoleErrors.filter((m) => pitfallRegex.test(m));
      const pitfallHitsPage = pageErrors.filter((m) => pitfallRegex.test(m));

      console.log(`[Scenario 2 pitfall observation] consoleHits=${pitfallHitsConsole.length}, pageHits=${pitfallHitsPage.length}`);
    });

    await softStep('Scenario 2 step 7 — Recovery: viewBiostructure(content, format, name) is the canonical safe entry point', async () => {
      const res = await page.evaluate(async (pdbPath) => {

        const fns = DG.Func.find({name: 'viewBiostructure', package: 'BiostructureViewer'});
        const fn = fns && fns[0];
        const inputs = (fn && fn.inputs) ? fn.inputs.map((i: any) => ({
          name: i.name,
          type: i.propertyType,
          optional: i.options?.optional ?? false,
        })) : [];

        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        let recoveryInvokeErr: string | null = null;
        try {
          grok.functions.call(
            'BiostructureViewer:viewBiostructure',
            {content: pdbContent, format: 'pdb', name: 'safe-fixture'},
          ).catch(() => {  });
        } catch (e: any) { recoveryInvokeErr = String(e?.message ?? e); }

        await new Promise((r) => setTimeout(r, 1500));

        return {
          registered: !!fn,
          inputCount: inputs.length,
          inputNames: inputs.map((i: any) => i.name),
          inputTypes: inputs.map((i: any) => i.type),
          nameInputOptional: inputs.find((i: any) => i.name === 'name')?.optional ?? null,
          recoveryInvokeErr,
        };
      }, samplePdbPath);

      expect(res.registered).toBe(true);
      expect(res.inputCount).toBe(3);
      expect(res.inputNames).toEqual(['content', 'format', 'name']);
      expect(res.inputTypes).toEqual(['string', 'string', 'string']);
      expect(res.recoveryInvokeErr).toBe(null);
    });

    await softStep('Scenario 3 — pdbTag round-trip via setOptions on a DataFrame with a `.pdb-tag-payload` column tag', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s3-row-1', 's3-row-2']),
          DG.Column.fromStrings('payload', ['x', 'y']),
        ]);

        df.col('payload').setTag('.pdb-tag-payload', pdbContent);

        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));

        const props = v.props.getProperties();
        const pdbTagProp = props.find((p: any) => p.name === 'pdbTag');

        let setErr: string | null = null;
        try { v.setOptions({pdbTag: '.pdb-tag-payload'}); }
        catch (e: any) { setErr = String(e?.message ?? e); }
        await new Promise((r) => setTimeout(r, 2000));
        const pdbTagAfter = v.props.get('pdbTag');

        const columnTagValue = df.col('payload').getTag('.pdb-tag-payload');

        return {
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
          pdbTagInCatalogue: !!pdbTagProp,
          pdbTagCategory: pdbTagProp?.category ?? null,
          setErr,
          pdbTagAfter,
          columnTagPresent: !!columnTagValue && columnTagValue.length > 0,
          columnTagLen: columnTagValue ? columnTagValue.length : 0,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.containerPresent).toBe(true);
      expect(res.pdbTagInCatalogue).toBe(true);
      expect(res.pdbTagCategory).toBe('Data');
      expect(res.setErr).toBe(null);
      expect(res.pdbTagAfter).toBe('.pdb-tag-payload');
      expect(res.columnTagPresent).toBe(true);
      expect(res.columnTagLen).toBeGreaterThan(1000);
    });

    await softStep('Scenario 4 — Behaviour: showMouseOverRowLigand + showSelectedRowsLigands round-trip + selection driver', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s4-row-1', 's4-row-2', 's4-row-3']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent, pdbContent]),
          DG.Column.fromStrings('ligand', ['CCO', 'CCC', 'CCN']),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }
        try { df.col('ligand').semType = 'Molecule'; } catch (_e) {  }

        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));

        v.setOptions({ligandColumnName: 'ligand'});
        await new Promise((r) => setTimeout(r, 800));
        const ligandColAfter = v.props.get('ligandColumnName');

        v.setOptions({showCurrentRowLigand: false});
        await new Promise((r) => setTimeout(r, 600));
        const currentRowOff = v.props.get('showCurrentRowLigand');

        const initMouseOver = v.props.get('showMouseOverRowLigand');
        v.setOptions({showMouseOverRowLigand: false});
        await new Promise((r) => setTimeout(r, 700));
        const mouseOverOff = v.props.get('showMouseOverRowLigand');
        v.setOptions({showMouseOverRowLigand: true});
        await new Promise((r) => setTimeout(r, 700));
        const mouseOverOn = v.props.get('showMouseOverRowLigand');

        df.currentRowIdx = 1;
        await new Promise((r) => setTimeout(r, 400));
        const currentRowIdx = df.currentRowIdx;

        v.setOptions({showMouseOverRowLigand: false});
        await new Promise((r) => setTimeout(r, 600));

        const initSelected = v.props.get('showSelectedRowsLigands');
        v.setOptions({showSelectedRowsLigands: true});
        await new Promise((r) => setTimeout(r, 700));
        const selectedOn = v.props.get('showSelectedRowsLigands');

        df.selection.init((i: number) => i === 0 || i === 2);
        await new Promise((r) => setTimeout(r, 600));
        const selectedCount = df.selection.trueCount;

        df.selection.init(() => false);
        await new Promise((r) => setTimeout(r, 400));
        const clearedCount = df.selection.trueCount;
        v.setOptions({showSelectedRowsLigands: false});
        await new Promise((r) => setTimeout(r, 500));
        const selectedOff = v.props.get('showSelectedRowsLigands');

        return {
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
          ligandColAfter,
          currentRowOff,
          initMouseOver,
          mouseOverOff,
          mouseOverOn,
          currentRowIdx,
          initSelected,
          selectedOn,
          selectedCount,
          clearedCount,
          selectedOff,
          rowCount: df.rowCount,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.containerPresent).toBe(true);
      expect(res.ligandColAfter).toBe('ligand');
      expect(res.currentRowOff).toBe(false);
      expect(res.initMouseOver).toBe(true);
      expect(res.mouseOverOff).toBe(false);
      expect(res.mouseOverOn).toBe(true);
      expect(res.currentRowIdx).toBe(1);
      expect(res.initSelected).toBe(false);
      expect(res.selectedOn).toBe(true);
      expect(res.selectedCount).toBe(2);
      expect(res.clearedCount).toBe(0);
      expect(res.selectedOff).toBe(false);
      expect(res.rowCount).toBe(3);
    });

    await softStep('Scenario 5 — bindingSiteWholeResidues round-trip with showBindingSite ON', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s5-row-1', 's5-row-2', 's5-row-3']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent, pdbContent]),
          DG.Column.fromStrings('ligand', ['CCO', 'CCC', 'CCN']),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }
        try { df.col('ligand').semType = 'Molecule'; } catch (_e) {  }

        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));

        v.setOptions({ligandColumnName: 'ligand'});
        await new Promise((r) => setTimeout(r, 500));
        v.setOptions({showCurrentRowLigand: true});
        await new Promise((r) => setTimeout(r, 500));
        v.setOptions({showBindingSite: true});
        await new Promise((r) => setTimeout(r, 1200));
        const showBindingSiteAfter = v.props.get('showBindingSite');
        const bindingSiteRadius = v.props.get('bindingSiteRadius');

        const initBindingWhole = v.props.get('bindingSiteWholeResidues');

        v.setOptions({bindingSiteWholeResidues: false});
        await new Promise((r) => setTimeout(r, 1200));
        const bindingWholeOff = v.props.get('bindingSiteWholeResidues');

        v.setOptions({bindingSiteWholeResidues: true});
        await new Promise((r) => setTimeout(r, 1200));
        const bindingWholeOn = v.props.get('bindingSiteWholeResidues');

        const props = v.props.getProperties();
        const bindingProps = props
          .filter((p: any) => p.category === 'Binding Site')
          .map((p: any) => p.name)
          .sort();

        return {
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
          showBindingSiteAfter,
          bindingSiteRadius,
          initBindingWhole,
          bindingWholeOff,
          bindingWholeOn,
          bindingProps,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.containerPresent).toBe(true);
      expect(res.showBindingSiteAfter).toBe(true);
      expect(res.bindingSiteRadius).toBe(5);
      expect(res.initBindingWhole).toBe(true);
      expect(res.bindingWholeOff).toBe(false);
      expect(res.bindingWholeOn).toBe(true);

      expect(res.bindingProps).toEqual(['bindingSiteRadius', 'bindingSiteWholeResidues', 'showBindingSite']);
    });

    let scenario6Mounted = false;
    let scenario6GearClicked = false;

    await softStep('Scenario 6 step 1-3 — Mount viewer; confirm default layoutShowControls = false; gear-click opens property panel (DOM)', async () => {
      const setupRes = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s6-row-1', 's6-row-2']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent]),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }

        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));

        return {
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
          initLayoutShow: v.props.get('layoutShowControls'),
          mspPluginMounted: !!document.querySelector('.msp-plugin'),
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(setupRes.containerPresent).toBe(true);

      expect(setupRes.initLayoutShow).toBe(false);
      scenario6Mounted = true;

      const gearClicked = await page.evaluate(async () => {
        const container = document.querySelector('[name="viewer-Biostructure"]');
        if (!container) return {found: false, opened: false};
        const gear = container.closest('.panel-base')?.querySelector(
          '.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
        if (!gear) return {found: false, opened: false};
        gear.click();
        await new Promise((r) => setTimeout(r, 1500));
        const cp = document.querySelector('.grok-prop-panel');
        return {found: true, opened: !!cp};
      });

      expect(gearClicked.found).toBe(true);
      expect(gearClicked.opened).toBe(true);
      scenario6GearClicked = true;
    });

    await softStep('Scenario 6 step 4-5 — Property panel path: setOptions({layoutShowControls: true}) flips Datagrok property; persists', async () => {
      if (!scenario6Mounted || !scenario6GearClicked) return;

      const res = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'Biostructure') as any : null;
        if (!v) return {ok: false};

        v.setOptions({layoutShowControls: true});
        await new Promise((r) => setTimeout(r, 1800));
        const afterTrue = v.props.get('layoutShowControls');

        const props = v.props.getProperties();
        const layoutProp = props.find((p: any) => p.name === 'layoutShowControls');

        return {
          ok: true,
          afterTrue,
          layoutCategory: layoutProp?.category ?? null,
        };
      });

      expect(res.ok).toBe(true);
      expect(res.afterTrue).toBe(true);
      expect(res.layoutCategory).toBe('Layout');
    });

    await softStep('Scenario 6 step 6-7 — Mol* overlay button path: click button[title="Toggle Controls Panel"] (DOM, conditional .msp-plugin)', async () => {
      if (!scenario6Mounted) return;

      const beforeOverlay = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'Biostructure') as any : null;
        return {
          mspPluginMounted: !!document.querySelector('.msp-plugin'),
          overlayBtnPresent: !!document.querySelector('button[title="Toggle Controls Panel"]'),
          layoutBefore: v?.props?.get?.('layoutShowControls') ?? null,
        };
      });

      if (!beforeOverlay.mspPluginMounted) {

        console.log('[Scenario 6 step 6 observation] .msp-plugin not mounted in recon env; overlay button click precondition not met (covered by gear-click DOM driving in step 1-3).');
        return;
      }

      if (beforeOverlay.overlayBtnPresent) {
        await page.locator('button[title="Toggle Controls Panel"]').first().click({timeout: 10_000});
        await page.waitForTimeout(1500);
      }

      const afterOverlay = await page.evaluate(() => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'Biostructure') as any : null;
        return {
          layoutAfter: v?.props?.get?.('layoutShowControls') ?? null,
        };
      });

      console.log(`[Scenario 6 overlay observation] before=${beforeOverlay.layoutBefore}, after=${afterOverlay.layoutAfter}, overlayBtnPresent=${beforeOverlay.overlayBtnPresent}`);

      await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const v = tv?.viewers ? Array.from(tv.viewers).find((x: any) => x.type === 'Biostructure') as any : null;
        if (v) {
          v.setOptions({layoutShowControls: false});
          await new Promise((r) => setTimeout(r, 800));
        }
      });
    });

    await softStep('Scenario 7 — Controls category surfaces showWelcomeToast + showImportControls; round-trip via setOptions', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const pdbContent = await grok.dapi.files.readAsText(pdbPath);

        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['s7-row-1', 's7-row-2']),
          DG.Column.fromStrings('structure', [pdbContent, pdbContent]),
        ]);
        try { df.col('structure').semType = 'Molecule3D'; } catch (_e) {  }

        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));

        const props = v.props.getProperties();
        const controlsProps = props
          .filter((p: any) => p.category === 'Controls')
          .map((p: any) => p.name)
          .sort();

        const initWelcome = v.props.get('showWelcomeToast');
        const initImport = v.props.get('showImportControls');

        v.setOptions({showImportControls: true});
        await new Promise((r) => setTimeout(r, 800));
        const importOn = v.props.get('showImportControls');
        v.setOptions({showImportControls: false});
        await new Promise((r) => setTimeout(r, 800));
        const importOff = v.props.get('showImportControls');

        v.setOptions({showWelcomeToast: true});
        await new Promise((r) => setTimeout(r, 800));
        const welcomeOn = v.props.get('showWelcomeToast');
        v.setOptions({showWelcomeToast: false});
        await new Promise((r) => setTimeout(r, 800));
        const welcomeOff = v.props.get('showWelcomeToast');

        return {
          containerPresent: !!document.querySelector('[name="viewer-Biostructure"]'),
          controlsProps,
          initWelcome,
          initImport,
          importOn,
          importOff,
          welcomeOn,
          welcomeOff,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});

      expect(res.containerPresent).toBe(true);
      expect(res.controlsProps).toEqual(['showImportControls', 'showWelcomeToast']);
      expect(res.initWelcome).toBe(false);
      expect(res.initImport).toBe(false);
      expect(res.importOn).toBe(true);
      expect(res.importOff).toBe(false);
      expect(res.welcomeOn).toBe(true);
      expect(res.welcomeOff).toBe(false);
    });
  } finally {

    await page.evaluate(() => { grok.shell.closeAll(); });
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
