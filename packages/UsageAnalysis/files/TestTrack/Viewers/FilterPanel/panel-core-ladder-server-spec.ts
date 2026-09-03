/* ---
realizes: [filters.cp.panel-core-ladder]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';
import {cardCount, clickResetCriteriaIcon, expectHeaderCounter, trueCount} from '../../helpers/filter-panel';
import {cardCaptions, checkboxCensus, driveHeaderSearch, establishTwoFilterState, filterState, FULL,
  headerSearchState, RACE_CATEGORY, visibleCardCaptions, waitForPanelSettled} from './panel-core-ladder-shared';

declare const grok: any;
declare const DG: any;

// Steps 10, 11, 11-negative and 13 of the ladder: the layout saved through the real View > Layout >
// Save to Gallery leaf and re-applied over the API, and the project round-trip. The gestures that
// build the two-filter state are proven in panel-core-ladder-spec.ts; here it is set through the
// filter group's API.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const projectName = `filters-ladder-${Date.now()}`;
const runTag = `ladder-${Date.now()}-${Math.floor(Math.random() * 1e9)}`;

async function myApplicableLayouts(page: Page): Promise<Array<{id: string; name: string}>> {
  return page.evaluate(async () => {
    const me = String(grok.shell.user.id);
    const ls = (await grok.dapi.layouts.getApplicable(grok.shell.tv.dataFrame)) ?? [];
    return ls
      .filter((l: any) => !l.author || !l.author.id || String(l.author.id) === me)
      .map((l: any) => ({id: String(l.id), name: String(l.friendlyName ?? l.name ?? '')}));
  });
}

function carriesMarker(name: string, marker: string): boolean {
  return name === marker || name.startsWith(`${marker} (`);
}

async function stampTableName(page: Page, marker: string): Promise<void> {
  const applied = await page.evaluate((name: string) => {
    grok.shell.tv.dataFrame.name = name;
    return String(grok.shell.tv.dataFrame.name);
  }, marker);
  expect(applied, 'the open table did not take this run\'s marker name, so a layout saved from it ' +
    'could not be told apart from one another process saved as this same user').toBe(marker);
}

// The marker is unique to this run, so the layouts applicable before the save need no listing.
async function saveLayoutToGallery(page: Page, marker: string): Promise<string> {
  const stamped = await page.evaluate(() => String(grok.shell.tv?.dataFrame?.name ?? ''));
  expect(stamped, 'the table is not carrying this run\'s marker at save time, so the layout the ' +
    'save produces would be unattributable').toBe(marker);

  expect(await v.driveTopMenuLeaf(page, ['View', 'Layout', 'Save to Gallery']),
    'the View | Layout | Save to Gallery leaf could not be driven, so no layout was saved').toBe(true);

  let mine: Array<{id: string; name: string}> = [];
  await expect.poll(async () => {
    mine = (await myApplicableLayouts(page)).filter((l) => carriesMarker(l.name, marker));
    return mine.length;
  }, {
    timeout: 25_000,
    intervals: [500, 1000, 2000, 3000],
    message: `the gallery save produced no new applicable layout named after this run's marker "${marker}"`,
  }).toBeGreaterThanOrEqual(1);
  if (mine.length !== 1) {
    throw new Error(
      `The gallery save produced ${mine.length} new applicable layouts carrying this run's marker ` +
      `"${marker}" (${mine.map((l) => `${l.id}:${l.name}`).join(', ')}), expected exactly 1 — ` +
      'refusing to guess which one is ours; deleting none.');
  }
  return mine[0].id;
}

async function deleteLayout(page: Page, layoutId: string | null): Promise<void> {
  if (!layoutId) return;
  await page.evaluate(async (id: string) => {
    try { const s = await grok.dapi.layouts.find(id); await grok.dapi.layouts.delete(s); } catch (_) {}
  }, layoutId);
}

test('Filter Panel — Panel Core Ladder: layout and project round-trips', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, withFilterPanel: true});
  await v.resetFilters(page);
  await expect.poll(() => trueCount(page), {timeout: 10_000}).toBe(FULL);
  await page.evaluate(() => grok.shell.tv.getFiltersGroup()
    .updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX'}));
  await expect.poll(() => cardCaptions(page), {timeout: 10_000}).toEqual(['SEX']);

  let trueCountSaved = -1;
  let projectId = '';

  try {
    await softStep('Step 10 Re-establish the two-filter state and save the layout, Step 11 layout ' +
      'round-trip → filters restored to the saved value', async () => {
      const mainMarker = `${runTag}-main`;
      let savedLayout: string | null = null;
      try {
        await stampTableName(page, mainMarker);
        trueCountSaved = await establishTwoFilterState(page);
        expect(trueCountSaved).toBeGreaterThan(0);
        expect(trueCountSaved).toBeLessThan(FULL);
        await expectHeaderCounter(page, '2', 'the two-filter state is re-established, so the counter must read 2 before the layout is saved');

        savedLayout = await saveLayoutToGallery(page, mainMarker);

        await page.evaluate(async () => {
          const was = grok.shell.tv.dataFrame.filter.trueCount;
          await (window as any).__filtered(() => {
            grok.shell.tv.dataFrame.resetFilter();
            try { grok.shell.tv.addViewer('Bar chart'); } catch (_) {}
          }, 800, was);
        });
        const perturbed = await trueCount(page);
        expect(perturbed).not.toBe(trueCountSaved);

        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          grok.shell.tv.loadLayout(saved);
        }, savedLayout);
        const settled = await waitForPanelSettled(page, {changedFrom: perturbed});
        const caps = await cardCaptions(page);
        expect(await page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]'))).toBe(true);
        expect(caps).toContain('RACE');
        expect(caps).toContain('AGE');
        const census = await checkboxCensus(page);
        expect(census.boxes).toBeGreaterThan(0);
        expect(census.boxes).toBe(census.cards);
        expect(census.checked).toBe(census.boxes);
        const raceRestored = await filterState(page, 'RACE', 'categorical');
        expect(raceRestored, 'the re-applied layout left the RACE card with no filter state')
          .not.toBeNull();
        expect(raceRestored!.selected,
          'the re-applied RACE card carries no category selection — the criterion did not survive')
          .not.toBeNull();
        expect(raceRestored!.selected).toEqual([RACE_CATEGORY]);
        const ageRestored = await filterState(page, 'AGE', 'histogram');
        expect(ageRestored, 'the re-applied layout left the AGE card with no filter state')
          .not.toBeNull();
        expect(ageRestored!.max, 'the re-applied AGE card carries no upper bound — its window was lost')
          .not.toBeNull();
        const ageColMaxAfterLayout = await page.evaluate(() => grok.shell.tv.dataFrame.col('AGE').max);
        expect(ageRestored!.max!).toBeLessThan(ageColMaxAfterLayout);
        expect(settled).toBe(trueCountSaved);
      } finally {
        await deleteLayout(page, savedLayout);
      }
    });

    await softStep('Step 11-negative GROK-16677 reset-then-save layout → full count, panel search empty', async () => {
      const resetMarker = `${runTag}-reset`;
      let resetLayout: string | null = null;
      try {
        await stampTableName(page, resetMarker);
        const cardsBefore = await cardCount(page);
        const countBeforeReset = await trueCount(page);
        expect(countBeforeReset,
          'nothing is filtering on the way into this reset, so the settle barrier below could not ' +
          'tell a completed reset from a stale reading').toBeLessThan(FULL);
        await clickResetCriteriaIcon(page, {via: 'dom'});
        const resetCount = await waitForPanelSettled(page, {changedFrom: countBeforeReset, timeoutMs: 10_000});
        expect(cardsBefore).toBeGreaterThan(0);
        expect(await cardCount(page)).toBe(cardsBefore);
        expect(resetCount).toBe(FULL);

        const searched = await driveHeaderSearch(page, 'RACE');
        expect(searched.typed).toBe('RACE');
        expect(searched.visibleAfter.length).toBeLessThan(searched.visibleBefore.length);
        expect(searched.countAfter).toBe(searched.countBefore);

        resetLayout = await saveLayoutToGallery(page, resetMarker);

        const perturbed = await page.evaluate(async () => {
          const was = grok.shell.tv.dataFrame.filter.trueCount;
          return (window as any).__filtered(() => grok.shell.tv.getFiltersGroup()
            .updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'RACE', selected: ['Asian']}), 900, was);
        });
        expect(perturbed).not.toBe(FULL);
        await page.evaluate(async (layoutId: string) => {
          const saved = await grok.dapi.layouts.find(layoutId);
          grok.shell.tv.loadLayout(saved);
        }, resetLayout);
        const settled = await waitForPanelSettled(page, {changedFrom: perturbed});
        expect(settled).toBe(FULL);
        const search = await headerSearchState(page);
        if (search === null) {
          throw new Error(
            'Step 11-negative: the re-applied panel has no header search input — its value cannot be read');
        }
        expect(search.value).toBe('');
        const census = await checkboxCensus(page);
        expect(census.cards).toBeGreaterThan(0);
        expect((await visibleCardCaptions(page)).length).toBe(census.cards);
      } finally {
        await deleteLayout(page, resetLayout);
      }
    });

    await softStep('Step 13 Project round-trip → panel reopens, filters restored (GROK-19152 barrier)', async () => {
      trueCountSaved = await establishTwoFilterState(page);

      const saved = await saveProjectViaApi(page, projectName);
      projectId = saved.projectId;

      await v.closeAllAndWait(page);

      await page.evaluate(async (id: string) => {
        const project = await grok.dapi.projects.find(id);
        await project.open();
      }, projectId);
      await page.waitForFunction(() => {
        const el = document.querySelector('[name="viewer-Filters"]') as HTMLElement | null;
        return !!el && el.offsetParent !== null;
      }, {timeout: 60_000});

      const reopened = await page.evaluate(async () => {
        const capsOf = () => Array.from(document.querySelectorAll('[name="viewer-Filters"] .d4-filter-column-name'))
          .map((c) => c.textContent?.trim());
        await (window as any).__poll(() => capsOf(),
          (c: string[]) => c.includes('RACE') && c.includes('AGE'), 1500, 50);
        const caps = capsOf();
        return {
          panel: !!document.querySelector('[name="viewer-Filters"]'),
          hasRace: caps.includes('RACE'),
          hasAge: caps.includes('AGE'),
          count: grok.shell.tv.dataFrame.filter.trueCount,
        };
      });
      expect(reopened.panel).toBe(true);
      expect(reopened.hasRace).toBe(true);
      expect(reopened.hasAge).toBe(true);
      const raceReopened = await filterState(page, 'RACE', 'categorical');
      expect(raceReopened, 'the reopened project left the RACE card with no filter state')
        .not.toBeNull();
      expect(raceReopened!.selected,
        'the reopened RACE card carries no category selection — the criterion was not persisted')
        .not.toBeNull();
      expect(raceReopened!.selected).toEqual([RACE_CATEGORY]);
      const ageReopened = await filterState(page, 'AGE', 'histogram');
      expect(ageReopened, 'the reopened project left the AGE card with no filter state').not.toBeNull();
      expect(ageReopened!.max, 'the reopened AGE card carries no upper bound — its window was not persisted')
        .not.toBeNull();
      const ageColMaxReopened = await page.evaluate(() => grok.shell.tv.dataFrame.col('AGE').max);
      expect(ageReopened!.max!).toBeLessThan(ageColMaxReopened);
      expect(reopened.count).toBe(trueCountSaved);
    });
  } finally {
    await deleteProjectWithCleanup(page, {projectId});
    await v.cleanupShell(page);
  }

  v.finishSpec();
});
