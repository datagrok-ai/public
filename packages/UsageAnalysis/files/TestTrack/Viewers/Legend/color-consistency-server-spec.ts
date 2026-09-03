/* ---
realizes: [viewers.histogram, viewers.line-chart]
--- */
import {expect, Page} from '@playwright/test';
import {test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {deleteEntities, layoutRoundTrip, projectRoundTrip} from './persistence';

// The server lane of the color-consistency scenario: the custom palette surviving a layout and a
// project round-trip. The grid coding, the picker and its propagation are in color-consistency-spec.ts
// on the local lane; the picked palette is set here directly.
test.use(specTestOptions);

const BLUE = '#1f77b4';
const probe = {column: 'Stereo Category', tags: ['.color-coding-categorical']};

function rOneLegendColors(page: Page): Promise<Record<string, string | null>> {
  return page.evaluate(() => {
    const out: Record<string, string | null> = {};
    for (const x of (window as any).grok.shell.tv.viewers) {
      if (x.type === 'Grid') continue;
      const items = Array.from(x.root.querySelectorAll('[name="legend"] .d4-legend-item')) as HTMLElement[];
      const rOne = items.find((el) => el.querySelector('.d4-legend-value')?.textContent?.trim() === 'R_ONE');
      out[x.type] = rOne ? getComputedStyle(rOne).color : null;
    }
    return out;
  });
}

const blueCount = (colors: Record<string, string | null>) =>
  Object.values(colors).filter((c) => c === 'rgb(31, 119, 180)').length;

test('Legend color consistency — palette persists across layout and project', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);
  await v.openTable(page);
  await v.addLegendViewers(page, {column: 'Stereo Category', viewers: ['Histogram', 'Line chart'], settleMs: 500});

  await softStep('Setup: custom palette R_ONE=blue, S_UNKN=green from the grid coding', async () => {
    const tag = await page.evaluate(async (blue) => {
      const col = (window as any).grok.shell.tv.dataFrame.col('Stereo Category');
      col.tags['.color-coding-type'] = 'Categorical';
      col.meta.colors.setCategorical({'R_ONE': blue, 'S_UNKN': '#00FF00'}, {fallbackColor: '#808080'});
      for (const x of (window as any).grok.shell.tv.viewers)
        if (x.type !== 'Grid') try { x.invalidate?.(); } catch (_) {}
      return String(JSON.parse(col.tags['.color-coding-categorical'] ?? '{}')['R_ONE'] ?? '').toLowerCase();
    }, BLUE);
    expect(tag).toBe(BLUE);
    expect(blueCount(await v.pollValue(() => rOneLegendColors(page), (c) => blueCount(c) >= 1, 3000, 100)),
      'at least 1 viewer renders R_ONE=blue before the round-trips').toBeGreaterThanOrEqual(1);
  });

  let layoutId: string | null = null;
  let projectId: string | null = null;
  try {
    await softStep('Save + re-apply layout — custom palette persists (tag verification)', async () => {
      const res = await layoutRoundTrip(page, 'ColorConsist', probe);
      layoutId = res.layoutId;
      expect(String(res.tags['.color-coding-categorical'] ?? '')).toBeTruthy();
      expect(String(JSON.parse(res.tags['.color-coding-categorical'] ?? '{}')['R_ONE'] ?? '').toLowerCase()).toBe(BLUE);
    });

    await softStep('Project round-trip — save + close + reopen + verify palette', async () => {
      const res = await projectRoundTrip(page, 'ColorConsistProj', probe);
      projectId = res.projectId;
      expect(res.ok, res.ok ? '' : `project save+reopen failed in phase '${res.phase}': ${res.error}`).toBe(true);
      expect(String(JSON.parse(res.tags['.color-coding-categorical'] ?? '{}')['R_ONE'] ?? '').toLowerCase()).toBe(BLUE);
      const colors = await v.pollValue(() => rOneLegendColors(page), (c) => blueCount(c) >= 1, 3000, 100);
      expect(blueCount(colors), 'at least 1 viewer reflects R_ONE=blue post-reopen').toBeGreaterThanOrEqual(1);
    });
  } finally {
    await softStep('Cleanup', async () => {
      await deleteEntities(page, {layoutIds: [layoutId], projectId});
    });
  }

  v.finishSpec();
});
