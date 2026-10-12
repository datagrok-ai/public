/* ---
sub_features_covered: [hmm.align.msa-engine, hmm.align.bio-dialog, hmm.homology.top-menu, hmm.logo.weblogo, hmm.build.top-menu]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';

test.use(specTestOptions);

const ANTIBODIES = 'System:AppData/Bio/samples/antibodies.csv';
const ROWS = 40;

/** Walks the pointer into a Bio top-menu group (d4 expands submenus off the mouse path). */
async function openBioMenu(page: Page, group: string, leaf: string): Promise<void> {
  const groupLoc = page.locator(`[name="div-Bio---${group}"]`);
  if (!await groupLoc.isVisible().catch(() => false))
    await page.locator('[name="div-Bio"]').click();
  await groupLoc.waitFor({state: 'visible', timeout: 15_000});
  const box = (await groupLoc.boundingBox())!;
  await page.mouse.move(box.x + 5, box.y + box.height / 2, {steps: 12});
  await page.mouse.move(box.x + box.width - 15, box.y + box.height / 2, {steps: 12});
  const leafLoc = page.locator(`[name="div-Bio---${group}---${leaf}"]`);
  await leafLoc.waitFor({state: 'visible', timeout: 15_000});
  await leafLoc.click();
}

test('HMM: align to profile in Bio\'s MSA dialog, homology search from the Bio menu', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async ({path, rows}) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const src = await grok.dapi.files.readCsv(path);
    const df = src.clone(DG.BitSet.create(src.rowCount, (i) => i < rows));
    df.name = 'antibodies';
    grok.shell.addTableView(df);
    for (let i = 0; i < 200; i++) {
      if (df.col('AntibodyHC')?.semType === 'Macromolecule' && document.querySelector('[name="viewer-Grid"] canvas'))
        break;
      await new Promise((r) => setTimeout(r, 200));
    }
  }, {path: ANTIBODIES, rows: ROWS});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 60_000});

  await softStep('S1: protein FASTA column offers Kalign and the hmmalign engine', async () => {
    await openBioMenu(page, 'Analyze', 'MSA...');
    const dlg = page.locator('[name="dialog-MSA"]');
    await dlg.waitFor({timeout: 60_000});
    const engine = dlg.locator('[name="input-host-Engine"] select');
    await engine.waitFor({state: 'visible', timeout: 30_000});
    const options = (await engine.locator('option').allTextContents()).filter((o) => o);
    expect(options[0]).toBe('Kalign');
    expect(options).toContain('Profile HMM (hmmalign)');
    await engine.selectOption({label: 'Profile HMM (hmmalign)'});
    const profile = dlg.locator('[name="input-host-Profile"] select');
    await profile.waitFor({state: 'visible', timeout: 30_000});
    await profile.selectOption({label: 'Antibody/TCR germline (auto)'});
    // The engine's long description wraps instead of widening the dialog.
    expect((await dlg.boundingBox())!.width).toBeLessThan(700);
    await page.screenshot({path: 'playwright/test-output/hmm-align-dialog.png'});
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.names().some((n) => n.startsWith('msa(')), null,
      {timeout: 120_000});
  });

  await softStep('S2: the aligned column keeps every residue, rows of equal length', async () => {
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const aligned = df.col(df.columns.names().find((n) => n.startsWith('msa('))!)!;
      const source = df.col('AntibodyHC')!;
      const lengths = new Set<number>();
      const changed: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        const row = aligned.get(i) ?? '';
        lengths.add(row.length);
        if (row.replace(/-/g, '') !== (source.get(i) ?? '').replace(/-/g, '').toUpperCase()) changed.push(i);
      }
      return {lengths: [...lengths], changed, units: aligned.meta.units, aligned: aligned.getTag('aligned')};
    });
    expect(info.lengths.length).toBe(1);
    expect(info.changed).toEqual([]);
    expect(info.units).toBe('fasta');
    expect(info.aligned).toContain('MSA');
  });
  await page.screenshot({path: 'playwright/test-output/hmm-align.png'});

  await softStep('S3: Homology Search (phmmer) with the current row as query', async () => {
    await page.evaluate(() => { grok.shell.tv.dataFrame.currentRowIdx = 2; });
    await openBioMenu(page, 'Search', 'Homology-Search-(phmmer)...');
    const dlg = page.locator('[name^="dialog-Homology-Search"]');
    await dlg.waitFor({timeout: 60_000});
    await dlg.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('phmmer row 3 E-value'), null,
      {timeout: 120_000});
    const info = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const e = df.col('phmmer row 3 E-value')!;
      const s = df.col('phmmer row 3 score')!;
      let best = -1;
      for (let i = 0; i < df.rowCount; i++) if (!s.isNone(i) && (best < 0 || s.get(i) > s.get(best))) best = i;
      return {reported: df.rowCount - e.stats.missingValueCount, best, selfE: e.get(2)};
    });
    expect(info.reported).toBeGreaterThan(1);
    expect(info.best).toBe(2); // the query's own row scores best
    expect(info.selfE).toBeLessThan(1e-50);
  });
  await page.screenshot({path: 'playwright/test-output/hmm-homology.png'});

  await softStep('S4: Profile HMM Logo adds a WebLogo of the chosen family', async () => {
    await openBioMenu(page, 'Analyze', 'Profile-HMM-Logo...');
    const dlg = page.locator('[name="dialog-Profile-HMM-Logo"]');
    await dlg.waitFor({timeout: 60_000});
    // The sequence column defaults to the hmmalign column of S1, whose germline model it proposes.
    await expect(dlg.locator('[name="input-host-Sequence"]')).toBeVisible();
    await dlg.locator('[name="input-host-Profile-from"] select').selectOption({label: 'Pfam family (built-in)'});
    await dlg.locator('[name="input-host-Family"] select').selectOption({label: 'C1-set'});
    await dlg.locator('[name="button-OK"]').click();
    await page.locator('[name="viewer-WebLogo"]').first().waitFor({timeout: 60_000});
    // Bio's WebLogo draws the profile HMM's positions (hmmlogo heights from this package).
    await page.waitForFunction(() => {
      const v = Array.from(grok.shell.tv.viewers).find((x) => x.type === 'WebLogo' &&
        (x as any).profile === 'C1-set') as any;
      return v?.profileData && v.positions?.length > 0;
    }, null, {timeout: 60_000});
    const positions = await page.evaluate(() => (Array.from(grok.shell.tv.viewers)
      .find((x) => x.type === 'WebLogo' && (x as any).profile === 'C1-set') as any).positions.length);
    expect(positions).toBe(86); // Pfam 37.0 C1-set has 86 match positions
    // Opened ready to read: a finite position width and a visible range (not a collapsed slider).
    await page.waitForTimeout(500);
    const layout = await page.evaluate(() => {
      const v = Array.from(grok.shell.tv.viewers).find((x) => x.type === 'WebLogo' &&
        (x as any).profile === 'C1-set') as any;
      return {width: v.positionWidth, shown: v.getWidgetStatus().values['positions shown']};
    });
    expect(Number.isFinite(layout.width)).toBe(true);
    expect(layout.shown).toBeGreaterThan(20);
  });
  await page.screenshot({path: 'playwright/test-output/hmm-logo.png'});

  await softStep('S5: Build Profile HMM from selected aligned rows, search the table, show its logo', async () => {
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i) => i < 15);
    });
    await openBioMenu(page, 'Search', 'Build-Profile-HMM...');
    const dlg = page.locator('[name="dialog-Build-Profile-HMM"]');
    await dlg.waitFor({timeout: 60_000});
    await dlg.locator('[name="input-host-Name"] input').fill('heavy15');
    await dlg.locator('[name="button-OK"]').click();
    // The new columns, or the error the dialog reports.
    await page.waitForFunction(() => grok.shell.tv.dataFrame.columns.contains('heavy15 E-value') ||
      Array.from(document.querySelectorAll('.d4-balloon-content')).some((b) => /failed/.test(b.textContent ?? '')),
    null, {timeout: 120_000});
    const error = await page.evaluate(() => Array.from(document.querySelectorAll('.d4-balloon-content'))
      .map((b) => b.textContent ?? '').find((t) => /failed/.test(t)));
    expect(error).toBeUndefined();
    const info = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const e = df.col('heavy15 E-value')!;
      let significant = 0;
      for (let i = 0; i < 15; i++) if (!e.isNone(i) && e.get(i) < 1e-10) significant++;
      for (let k = 0; k < 150; k++) {
        const v = Array.from(grok.shell.tv.viewers).find((x) => x.type === 'WebLogo' &&
          (x as any).profile === 'heavy15') as any;
        if (v?.profileData && v.positions?.length > 0) {
          return {significant, logo: v.positions.length, mapped: v.profileData.columnPositions?.length ?? 0,
            kept: !!df.getTag('.hmm.profile.heavy15')};
        }
        await new Promise((r) => setTimeout(r, 100));
      }
      return {significant, logo: 0, mapped: 0, kept: false};
    });
    expect(info.significant).toBe(15); // every training row is found
    expect(info.logo).toBeGreaterThan(100);
    expect(info.kept).toBe(true); // the profile is kept with the table
    expect(info.mapped).toBe(info.logo); // its positions map to the alignment it was built from

    // Clicking the tallest letter selects the alignment's rows that have it at that position.
    const target = await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.setAll(false);
      const v = Array.from(grok.shell.tv.viewers).find((x) => x.type === 'WebLogo' &&
        (x as any).profile === 'heavy15') as any;
      const status = v.getWidgetStatus();
      const areas = Object.entries(status.hitAreas as Record<string, {x: number; y: number; width: number;
        height: number}>).filter(([k]) => k.startsWith('monomer '));
      const [name, a] = areas.sort(([, p], [, q]) => q.height - p.height)[0];
      const r = (status.parts.canvas as HTMLCanvasElement).getBoundingClientRect();
      return {name, x: r.left + a.x + a.width / 2, y: r.top + a.y + a.height / 2};
    });
    await page.mouse.click(target.x, target.y);
    const selected = await page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount);
    expect(selected, `rows selected by ${target.name}`).toBeGreaterThan(0);
  });
  await page.screenshot({path: 'playwright/test-output/hmm-build.png'});
  finishSpec(stepErrors);
});
