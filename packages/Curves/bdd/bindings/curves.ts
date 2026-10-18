/* Curves' own scene setting. The datasets are the package's files on the stand (`System:AppData/Curves/...`,
   published with the package). The legend of the Multi curve viewer is read through the viewer's own
   status (`getWidgetStatus`: a `legend row N` / `shortened legend row N` area per row, `legend more`, and
   the `curves` and `legend …` readings), so the features use the library's area and reading steps. */
import {Page} from '@playwright/test';
import {dataset, Given} from '@datagrok-libraries/bdd';
import {silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

dataset('assay-curves', {path: 'System:AppData/Curves/demo/Assay_Curves.csv', aliases: ['Assay Curves'],
  description: '43 dose-response curves in "Fitted Curve", one per row, with "Compound|Assay|Target" naming each'});

/* The legend needs a plot of 325 x 275 px at least (`FitConstants.MIN_LEGEND_PX_*`): beside the grid in a
   help GIF's 800 px window the viewer docked at its default half would draw none, and no library step
   docks a viewer at a share of the view. It overlays the selected rows only, as in the Assay Curves demo:
   the current row's curve is one of them, and would be named twice. */
export const addWideViewer = Given('user adds a MultiCurveViewer viewer of the selected rows taking {int} percent of the view width', async (page: Page, percent: number) => {
  silent(page);
  await page.evaluate(async (ratio) => {
    const tv = grok.shell.tv;
    const v = await tv.dataFrame.plot.fromType('MultiCurveViewer', {showCurrentRowCurve: false, showMouseOverRowCurve: false});
    tv.dockManager.dock(v, DG.DOCK_TYPE.RIGHT, null, 'Multi curve viewer', ratio);
  }, percent / 100);
  await page.locator('[name="viewer-MultiCurveViewer"]').filter({visible: true}).first().waitFor();
}, {tier: 'api', description: 'a multi curve viewer on the current table, docked to the right of the grid at that share of the view; silent in a guide'});

/* A grid cell draws a legend only from 325 x 275 px inside its 6 px margin, wider than a grid makes a
   column by default; no library step sets a grid column's width. */
export const columnWidth = Given('the {string} column of the grid is {int} px wide', async (page: Page, column: string, width: number) => {
  silent(page);
  await page.evaluate(([name, w]) => { grok.shell.tv.grid.columns.byName(name).width = w; }, [column, width] as [string, number]);
}, {tier: 'api', description: 'the width of one column of the current table view\'s grid; silent in a guide'});
