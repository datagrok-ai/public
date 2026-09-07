/* ---
realizes: [filters.cp.save-and-reapply-state]
--- */
import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {dismissPanelMenu, readSaveOrApplyLeaves, removeProbeState, saveStateViaMenu} from './save-and-reapply-shared';

declare const grok: any;

// Step 8 of the scenario: the demog-typed state must not be offered on beer.csv, which has no
// local copy. The state is saved here the same way Step 4 saves it (the localStorage of this
// lane is not the local lane's), on a RACE / AGE configuration set up through the API.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const beerPath = 'System:DemoFiles/beer.csv';
const ALWAYS_OFFERED_LEAF = 'Save...';
const probeName = `filter-state-${Date.now()}`;

test('Filters — a saved filter state is not offered on a mismatched table', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true, semTypeTimeoutMs: 1000});

  try {
    await softStep('Setup — RACE = Asian and AGE 30-60 filtering, state saved through Save or Apply > Save...', async () => {
      await v.resetFilters(page);
      const {filteredCount} = await v.applyCategoricalFilter(page, 'RACE', ['Asian']);
      const count = await v.applyNumericFilter(page, 'AGE', 30, 60);
      expect(filteredCount).toBeGreaterThan(0);
      expect(count).toBeGreaterThan(0);
      expect(count).toBeLessThan(filteredCount);
      await saveStateViaMenu(page, probeName);
    });

    await softStep('Step 8 Probe name absent under Save or Apply on beer table', async () => {
      await v.openTable(page, {path: beerPath, withFilterPanel: true, semTypeTimeoutMs: 1000});
      await expect.poll(async () => page.evaluate(() =>
        (grok.shell.tv.root as HTMLElement)
          .querySelectorAll('[name="viewer-Filters"] .d4-filter-group-header').length), {
        timeout: 20_000,
        intervals: [30, 60, 120, 250, 500, 1000],
        message: 'the beer table view never grew a Filter Panel to read the Save or Apply submenu from',
      }).toBe(1);
      const leaves = await readSaveOrApplyLeaves(page);
      expect(leaves,
        'the Save or Apply submenu enumerated nothing on beer.csv, so the absence of the probe name ' +
        'below would be satisfied by a submenu that simply failed to populate')
        .toContain(ALWAYS_OFFERED_LEAF);
      expect(leaves,
        'the demog-typed state is offered on beer.csv, whose columns do not match it')
        .not.toContain(probeName);
      await dismissPanelMenu(page);
    });
  } finally {
    await removeProbeState(page, probeName);
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
