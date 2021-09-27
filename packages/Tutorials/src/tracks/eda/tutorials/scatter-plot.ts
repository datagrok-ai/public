import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class ScatterPlotTutorial extends Tutorial {
  get name() { return 'Scatter Plot'; }
  get description() {
    return 'A graph in which the values of two variables are plotted along two axes';
  }
  get steps(){ return 3} //set thew number of steps
    
  protected async _run() {
    this.header.textContent = 'Scatter Plot';
    const plot = <DG.ScatterPlotViewer>(await this.openPlot(
      'scatter plot',
      (x) => x.type === DG.VIEWER.SCATTER_PLOT));
    const info = <{ [key: string]: any }>plot.getInfo();
    const columnCheck = (selector: DG.ColumnComboBox, col: string) => selector.onChanged.pipe(filter((v) => v === col));

    this.describe("Now, let's customize it by choosing the columns of interest on corresponding axes. " +
      "The easiest way to do it is to click on the column selector on the viewer. Alternatively, " +
      "you can drag the column right from the spreadsheet, or from the column list " +
      "(Windows | Columns, or Alt+C). Please try different ways to do it.");

    await this.action('Set X to HEIGHT', columnCheck(info.xColSelector, 'HEIGHT'), info.xColSelector.root);
    await this.action('Set Y to WEIGHT', columnCheck(info.yColSelector, 'WEIGHT'), info.yColSelector.root);
    await this.action('Set Size to AGE', columnCheck(info.sizeColSelector, 'AGE'), info.sizeColSelector.root);
    await this.action('Set Color to SEX', columnCheck(info.colorColSelector, 'SEX'), info.colorColSelector.root);

    this.describe("To zoom in, hold Alt key, and drag a rectangle that you want to zoom in to. Try it now.");
    await this.action('Zoom in', plot.onZoomed);

    this.describe("As in most viewers, double-clicking on an empty area resets the view. Alternatively, " +
      "this option is always available in the context menu.");
    await this.action('Double-click to unzoom', plot.onResetView);

    this.describe("Click on a point. Note that it becomes the current point in a spreadsheet, too.");
    await this.action('Click on a point', this.t.onCurrentRowChanged);

    this.describe("Select points by dragging a rectangle on a viewer while holding Shift. " +
      "Note that the row selection is being reflected on most viewers, such as the spreadsheet. " +
      "A number of points under the selection rectangle is shown right there.");
    await this.action('Select points', this.t.onSelectionChanged);

    this.describe("Deselect some points by dragging a rectangle on a viewer while holding Ctrl+Shift.");
    await this.action('Deselect points', this.t.onSelectionChanged);
  }
}
