/* The one step the correlation plot needs beyond the library's `viewers` tier.

   Everything else its features say is read from what the plot reports about itself
   (`cell <X> x <Y>` / `row header <Y>` / `pinned band` regions, `correlation of <X> and <Y>`,
   `text of cell`, `color of cell`, `cell type of`, `last clicked cell` / `last clicked value`,
   `tooltip plot x column` — see `core/client/d4/lib/src/viewers/correlation_plot/CLAUDE.md`), so
   the probe-click calibration loop, the five-pixel median colour sampler and the nine
   `cp.getCorrelation(...)` reflections of the three TestTrack specs have no successor here.

   What is left is the one claim the plot cannot make about itself: that the number in a cell is
   the coefficient of THOSE TWO COLUMNS over THOSE ROWS, computed by something other than the
   viewer. `DG.Stats` is that something; the pivot table's `aggregationMatches` is the precedent. */
import {expect, Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';
import {ElementRef, viewers} from '@datagrok-libraries/bdd/runtime';

/** The coefficient one cell holds, against `DG.Stats` over the rows the viewer's Row Source names.
 * Both are pairwise-complete, so a pair with blanks is compared over the same rows on both sides;
 * the reference is built from the source table, never from the viewer. The step refuses a viewer
 * whose own `filter` formula is set — that would mean re-implementing the formula here, and the
 * feature states the measured value instead. */
export const correlationMatches = Then(
  'the correlation of {string} and {string} of {widget} should match the {word} coefficient of the table',
  async (page: Page, x: string, y: string, target: ElementRef, kind: string) => {
    await viewers.installViewerRuntime(page);
    const loc = await viewers.viewerLocator(page, target);
    const report: string = await loc.evaluate((element, [xName, yName, want]) => {
      const bdd = (window as any).__bdd;
      const DG = (window as any).DG;
      const viewer = bdd.viewerOf(element);
      const values = viewer.getWidgetStatus().values;
      const key = `correlation of ${xName} and ${yName}`;
      if (!(key in values)) {
        const held = Object.keys(values).filter((k) => k.startsWith('correlation of ')).join('; ');
        return `the matrix holds no "${key}"; it holds: ${held}`;
      }
      const type = String(values['correlation type'] ?? '');
      if (type.toLowerCase() !== want.toLowerCase())
        return `the matrix is showing ${type} coefficients, not ${want}`;
      const formula = String(viewer.props.filter ?? '');
      if (formula !== '')
        return `the viewer's own filter is "${formula}" — this check reads the table's rows, not a formula`;
      const source = viewer.dataFrame;
      const rowSource = String(viewer.props.rowSource ?? 'Filtered');
      const rows = rowSource === 'Filtered' ? source.filter
        : rowSource === 'Selected' ? source.selection
        : rowSource === 'All' ? null : undefined;
      if (rows === undefined)
        return `Row Source is "${rowSource}" — this check knows All, Filtered and Selected`;
      const over = rows === null ? source : source.clone(rows);
      const stats = DG.Stats.fromColumn(over.col(xName));
      const reference = want.toLowerCase() === 'spearman'
        ? stats.spearmanCorr(over.col(yName)) : stats.corr(over.col(yName));
      const shown = Number(values[key]);
      if (!isFinite(shown) || !isFinite(reference))
        return `the matrix shows ${values[key]} and ${want} over ${over.rowCount} rows is ${reference}`;
      if (Math.abs(shown - reference) > 1e-6) {
        return `${xName}/${yName}: the matrix shows ${shown}, ${want} over the ${over.rowCount} `
          + `${rowSource.toLowerCase()} rows is ${reference} (pearson source: ${values['pearson source']})`;
      }
      return '';
    }, [x, y, kind] as [string, string, string]);
    expect(report, `${target.phrase}`).toBe('');
  }, {description: 'the number in the cell against DG.Stats over the rows the Row Source names — the only check that does not ask the viewer twice'});
