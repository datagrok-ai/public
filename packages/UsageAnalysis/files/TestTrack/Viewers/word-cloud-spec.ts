import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Word cloud', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    const df = await w.grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
    w.grok.shell.addTableView(df);
    await new Promise((resolve: any) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const cols: any[] = [];
    for (let i = 0; i < df.columns.length; i++) cols.push(df.columns.byIndex(i));
    const hasBioChem = cols.some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r: any) => setTimeout(r, 200));
      }
      await new Promise((r: any) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Step 2: Open Add Viewer gallery, click Word Cloud tile', async () => {
    await page.locator('i[aria-label="Add viewer"]').click();
    await page.waitForTimeout(800);
    const added = await page.evaluate(async () => {
      const w = window as any;
      const cards = Array.from(document.querySelectorAll('.d4-item-card.viewer-gallery'))
        .filter(e => !!(e as HTMLElement).offsetParent);
      const target = cards.find(e => (e as HTMLElement).innerText.trim() === 'Word Cloud') as HTMLElement | undefined;
      if (!target) return false;
      ['pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click'].forEach(type => {
        target.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, view: window, button: 0}));
      });
      await new Promise((r: any) => setTimeout(r, 1500));
      return !!w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(added).toBe(true);
  });

  await softStep('Step 3: Close viewer, re-add via Toolbox Word Cloud icon', async () => {
    const reAdded = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Word-cloud"]');
      if (root) {
        const panel = root.closest('.panel-base') || root.parentElement!;
        const closeBtn = panel.querySelector('[name="Close"]') || panel.querySelector('[name="icon-times"]');
        (closeBtn as HTMLElement)?.click();
        await new Promise((r: any) => setTimeout(r, 500));
      }
      const toolboxIcon = document.querySelector('[name="icon-Word-cloud"]') as HTMLElement | null;
      if (!toolboxIcon) return false;
      ['pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click'].forEach(type => {
        toolboxIcon.dispatchEvent(new MouseEvent(type, {bubbles: true, cancelable: true, view: window, button: 0}));
      });
      await new Promise((r: any) => setTimeout(r, 1500));
      return !!w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(reAdded).toBe(true);
  });

  await softStep('Step 4: Hamburger menu opens (canvas text clicks not DOM-testable)', async () => {
    const menuText = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      const panel = root.closest('.panel-base') || root.parentElement!;
      const hamburger = panel.querySelector('[name="icon-font-icon-menu"]') as HTMLElement;
      hamburger.click();
      await new Promise((r: any) => setTimeout(r, 600));
      const popup = Array.from(document.querySelectorAll('.d4-menu-popup'))
        .find(e => !!(e as HTMLElement).offsetParent) as HTMLElement | undefined;
      const text = popup ? popup.innerText : '';
      document.body.click();
      return text;
    });
    expect(menuText).toContain('Properties');
  });

  await softStep('Step 5: Open Property Pane via Gear icon', async () => {
    const hasData = await page.evaluate(async () => {
      const root = document.querySelector('[name="viewer-Word-cloud"]')!;
      const panel = root.closest('.panel-base') || root.parentElement!;
      const gear = panel.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      gear.click();
      await new Promise((r: any) => setTimeout(r, 800));
      return !!Array.from(document.querySelectorAll('.property-grid-category'))
        .find(e => e.textContent?.trim() === 'Data');
    });
    expect(hasData).toBe(true);
  });

  await softStep('Step 6: Modify column, minTextSize, maxTextSize, bold, rotationStep', async () => {
    const after = await page.evaluate(async () => {
      const w = window as any;
      const wc = w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
      wc.props.columnColumnName = 'Stereo Category';
      await new Promise((r: any) => setTimeout(r, 400));
      wc.props.minTextSize = 20;
      wc.props.maxTextSize = 80;
      wc.props.bold = true;
      wc.props.rotationStep = 90;
      await new Promise((r: any) => setTimeout(r, 600));
      return {
        column: wc.props.columnColumnName,
        minTextSize: wc.props.minTextSize,
        maxTextSize: wc.props.maxTextSize,
        bold: wc.props.bold,
        rotationStep: wc.props.rotationStep,
      };
    });
    expect(after.column).toBe('Stereo Category');
    expect(after.minTextSize).toBe(20);
    expect(after.maxTextSize).toBe(80);
    expect(after.bold).toBe(true);
    expect(after.rotationStep).toBe(90);
  });

  await softStep('Cleanup: Close viewer', async () => {
    const gone = await page.evaluate(async () => {
      const w = window as any;
      const root = document.querySelector('[name="viewer-Word-cloud"]');
      if (root) {
        const panel = root.closest('.panel-base') || root.parentElement!;
        const closeBtn = panel.querySelector('[name="Close"]') || panel.querySelector('[name="icon-times"]');
        (closeBtn as HTMLElement)?.click();
        await new Promise((r: any) => setTimeout(r, 500));
      }
      return !w.grok.shell.tv.viewers.find((v: any) => v.type === 'Word cloud');
    });
    expect(gone).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
});
