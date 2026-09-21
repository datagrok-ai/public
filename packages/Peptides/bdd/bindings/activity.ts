import {expect, type Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {el, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

function expectScaled(source: number[], actual: number[], scaling: string): void {
  expect(['none', 'lg', '-lg']).toContain(scaling);
  expect(source.length, 'source activities are populated').toBeGreaterThan(0);
  expect(actual.length).toBe(source.length);
  source.forEach((value, row) => {
    expect(Number.isFinite(value) && value > 0, `source activity of row ${row + 1}`).toBe(true);
    const expected = scaling === 'none' ? value : Math.log10(value) * (scaling === '-lg' ? -1 : 1);
    // Activity is a Float32 column; allow its rounding, not a different transformation.
    expect(Math.abs(actual[row] - expected), `scaled activity of row ${row + 1}`)
      .toBeLessThanOrEqual(Math.max(Math.abs(expected), 1) * 1e-7);
  });
}

export const previewScaling = Then('the peptide activity preview should use {string} scaling',
  async (page: Page, scaling: string) => {
    const target = el('Histogram viewer in Peptides pane');
    await viewers.settle(page, target);
    const values = await viewers.onViewer(page, target, (root) => {
      const histogram = DG.Widget.find(root);
      return {source: grok.shell.t.getCol('IC50').toList(), actual: histogram.dataFrame.getCol('Activity').toList()};
    });
    expectScaled(values.source, values.actual, scaling);
  }, {description: 'compares every plotted Activity value with an independent transform of source IC50'});

export const analysisScaling = Then('the SAR activity column should use {string} scaling',
  async (page: Page, scaling: string) => {
    const values = await page.evaluate(() => ({
      source: grok.shell.t.getCol('IC50').toList(), actual: grok.shell.t.getCol('Activity').toList(),
    }));
    expectScaled(values.source, values.actual, scaling);
  }, {description: 'compares the actual analysis column against every source IC50 value'});
