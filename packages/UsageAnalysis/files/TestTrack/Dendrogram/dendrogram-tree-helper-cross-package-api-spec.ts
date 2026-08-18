import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Dendrogram: Cross-package TreeHelper / DendrogramService consumption (JS API smoke)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await softStep('Scenario 1.1 - Dendrogram:getTreeHelper resolves to an ITreeHelper with the documented method surface', async () => {
    const result = await page.evaluate(async () => {
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');
      const required = ['newickToDf', 'toNewick', 'getLeafList', 'getNodeList'];
      const methods = required.map(m => ({m, isFn: typeof th?.[m] === 'function'}));
      return {
        truthy: th != null,
        constructorName: th?.constructor?.name,
        methods,
        allMethodsCallable: methods.every(x => x.isFn),
      };
    });
    expect(result.truthy, 'getTreeHelper returns a non-null ITreeHelper').toBe(true);
    expect(result.allMethodsCallable,
      `getTreeHelper instance exposes the ITreeHelper method surface (${result.methods.map(x => `${x.m}=${x.isFn}`).join(', ')})`)
      .toBe(true);
  });

  await softStep('Scenario 1.2 - Dendrogram:getDendrogramService resolves to an IDendrogramService with callable injectTreeForGrid', async () => {
    const result = await page.evaluate(async () => {
      const svc: any = await grok.functions.call('Dendrogram:getDendrogramService');
      return {
        truthy: svc != null,
        constructorName: svc?.constructor?.name,
        hasInject: typeof svc?.injectTreeForGrid === 'function',
      };
    });
    expect(result.truthy, 'getDendrogramService returns a non-null IDendrogramService').toBe(true);
    expect(result.hasInject, 'IDendrogramService.injectTreeForGrid is a callable function').toBe(true);
  });

  await softStep('Scenario 1.3 - injectTreeForGrid(grid, 4-leaf treeRoot, "leaf") mounts the grid neighbor with the magic-wand icon and emits no fatal console error', async () => {
    const result = await page.evaluate(async () => {

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 600));

      const df = DG.DataFrame.fromColumns([
        DG.Column.fromList('string', 'leaf', ['A', 'B', 'C', 'D']),
      ]);
      df.name = 'tree-helper-cross-package-smoke';

      const tv = grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1200));

      const treeRoot: any = {
        name: 'root',
        branch_length: 0,
        children: [
          {name: 'I1', branch_length: 1, children: [
            {name: 'A', branch_length: 1, children: []},
            {name: 'B', branch_length: 1, children: []},
          ]},
          {name: 'I2', branch_length: 1, children: [
            {name: 'C', branch_length: 1, children: []},
            {name: 'D', branch_length: 1, children: []},
          ]},
        ],
      };

      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');

      const svc: any = await grok.functions.call('Dendrogram:getDendrogramService');

      const consoleErrors: string[] = [];
      const origErr = console.error;
      console.error = (...args: any[]) => {
        consoleErrors.push(args.map((a: any) => String(a)).join(' '));
        origErr.apply(console, args);
      };

      let injectThrew: string | false = false;
      try {
        svc.injectTreeForGrid(tv.grid, treeRoot, 'leaf');
      } catch (e: any) {
        injectThrew = String(e?.message || e);
      }

      let magicWandMounted = false;
      let closeBttnMounted = false;
      for (let i = 0; i < 30; i++) {
        if (document.querySelector('.dendrogram-assign-clusters-bttn')) {
          magicWandMounted = true;
          closeBttnMounted = !!document.querySelector('.dendrogram-close-bttn');
          break;
        }
        await new Promise(r => setTimeout(r, 300));
      }

      console.error = origErr;

      const fatalErrors = consoleErrors.filter(t =>
        !/Failed to load resource[\s\S]*404/i.test(t)
        && !/ResizeObserver loop/i.test(t));

      const viewerTypes = Array.from(grok.shell.tv?.viewers || []).map((v: any) => v.type);

      return {
        injectThrew,
        magicWandMounted,
        closeBttnMounted,
        fatalErrors,
        viewerTypes,
      };
    });

    expect(result.injectThrew,
      'svc.injectTreeForGrid(grid, treeRoot, "leaf") must NOT throw')
      .toBe(false);
    expect(result.magicWandMounted,
      '.dendrogram-assign-clusters-bttn (magic-wand) mounts on the grid neighbor within budget - confirms injectTreeForGrid attached the GridNeighbor (atlas: dendrogram.api.dendrogram-service.inject-tree-for-grid)')
      .toBe(true);
    expect(result.closeBttnMounted,
      '.dendrogram-close-bttn also mounts (companion neighbor affordance documented in grok-browser/references/dendrogram.md)')
      .toBe(true);
    expect(result.fatalErrors,
      'no fatal console errors emitted during getTreeHelper + getDendrogramService + injectTreeForGrid sequence')
      .toEqual([]);
    expect(result.viewerTypes,
      'TableView viewers list still enumerates [Grid] only (the dendrogram neighbor is a GridNeighbor, not a DG.Viewer)')
      .toEqual(['Grid']);
  });

  await softStep('Scenario 1 -> Scenario 2 transition: close the dendrogram neighbor from Scenario 1.3', async () => {
    await page.evaluate(async () => {
      const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
      if (closeBtn) {
        closeBtn.click();
        await new Promise(r => setTimeout(r, 400));
      }
      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));
    });
  });

  await softStep('Scenario 2 - getDendrogramService singleton invariant: svc2 === svc1, window.$dendrogramService === svc1 on first AND second call', async () => {
    const result = await page.evaluate(async () => {

      const initialCached: any = (window as any).$dendrogramService;
      const initialWasSet = initialCached != null;

      const svc1: any = await grok.functions.call('Dendrogram:getDendrogramService');

      const winAfter1 = (window as any).$dendrogramService;
      const win1EqualsSvc1 = winAfter1 === svc1;
      const winAfter1Set = winAfter1 != null;

      const svc2: any = await grok.functions.call('Dendrogram:getDendrogramService');

      const svc2EqualsSvc1 = svc2 === svc1;

      const winAfter2 = (window as any).$dendrogramService;
      const win2EqualsSvc1 = winAfter2 === svc1;

      const svc1RespectsInitial = !initialWasSet || (svc1 === initialCached);

      return {
        initialWasSet,
        winAfter1Set,
        win1EqualsSvc1,
        svc2EqualsSvc1,
        win2EqualsSvc1,
        svc1RespectsInitial,
        svc1HasInject: typeof svc1?.injectTreeForGrid === 'function',
      };
    });
    expect(result.winAfter1Set,
      'window.$dendrogramService is set after the first getDendrogramService call (atlas: dendrogram.api.get-dendrogram-service singleton)')
      .toBe(true);
    expect(result.win1EqualsSvc1,
      'window.$dendrogramService === svc1 after the first call (referential equality)')
      .toBe(true);
    expect(result.svc2EqualsSvc1,
      'svc2 === svc1 (singleton instance is reused, not re-created on the second call)')
      .toBe(true);
    expect(result.win2EqualsSvc1,
      'window.$dendrogramService === svc1 continues to hold after the second call')
      .toBe(true);
    expect(result.svc1RespectsInitial,
      'when window.$dendrogramService pre-existed at Scenario 2 step 1, svc1 equals that cached instance (the package.ts#L84 guard returns the cache, does not overwrite)')
      .toBe(true);
    expect(result.svc1HasInject,
      'the singleton instance exposes the IDendrogramService surface (callable injectTreeForGrid)')
      .toBe(true);
  });

  await page.evaluate(() => {
    const closeBtn = document.querySelector('.dendrogram-close-bttn') as HTMLElement | null;
    if (closeBtn) closeBtn.click();
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
