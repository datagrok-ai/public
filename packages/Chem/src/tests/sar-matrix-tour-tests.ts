import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {category, delay, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {dockSarMatrixTabs} from '../analysis/sar-matrix/sar-matrix-viewer';
import {addSarMatrixTour, SarHint} from '../demo/sar-matrix-tour';


category('SAR Matrix: demo tour', () => {
  test('shown once, out of sight with its view, started over by the icon', async () => {
    const demoName = 'TourTest';
    const shownKey = `chem-demo-tour-shown-${demoName}`;
    localStorage.removeItem(shownKey);
    const tv = grok.shell.addTableView(grok.data.demo.demog(10));
    let other: DG.TableView | null = null;
    try {
      const viewer = tv.scatterPlot();
      dockSarMatrixTabs(tv, viewer);
      viewer.root.append(ui.div([], 'chem-sar-sum-overview')); // an analysis on screen is what a tour waits for
      const hints: SarHint[] = ['First', 'Second'].map((title) => ({anchor: () => viewer.root, title, text: title}));
      const cards = (): string =>
        [...document.querySelectorAll('.chem-sar-tour-card h3')].map((h) => h.textContent).join(', ');
      const shows = (expected: string): Promise<void> =>
        waitUntil(() => cards() === expected, () => `the cards are "${cards()}", not "${expected}"`);
      const click = (selector: string): void => document.querySelector<HTMLElement>(selector)!.click();

      addSarMatrixTour(viewer, demoName, hints);
      await shows('First');
      click('.chem-sar-tour-card button');
      await shows('Second');
      other = grok.shell.addTableView(grok.data.demo.demog(10));
      await shows('');
      grok.shell.v = tv;
      await shows('Second');
      click('.chem-sar-tour-card .fa-times');
      await shows('');

      addSarMatrixTour(viewer, demoName, hints);
      await delay(500);
      expect(cards(), '', 'the cards of a demo opened again');
      click('.chem-sar-tour-icon');
      await shows('First');
      tv.close();
      await shows('');
    } finally {
      localStorage.removeItem(shownKey);
      tv.close();
      other?.close();
    }
  });
});
