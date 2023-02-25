import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';


export class ScatterPlotTutorial extends Tutorial {
  get name() {
    return 'Scatter Plot';
  }
  get description() {
    return 'A graph in which the values of two variables are plotted along two axes';
  }
  get steps() { return 11; }
  
  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers/scatter-plot';

  protected async _run() {
    this.header.textContent = this.name;

    this.describe('A scatter plot is a mathematical diagram that uses Cartesian coordinates ' +
      'to display values for typically two variables for a set of data. If the points are ' +
      'encoded by size, color, or marker shape, you can increase the number of displayed variables. ' +
      'The data is displayed as a collection of points, each having the value of one variable ' +
      'determining the position on the horizontal axis and the value of the other variable ' +
      'determining the position on the vertical axis.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    const plot = (await this.openPlot(
      'scatter plot',
      (x) => x.type === DG.VIEWER.SCATTER_PLOT)) as unknown as DG.ScatterPlotViewer;
    const info = <{ [key: string]: any }>plot.getInfo();
    const columnCheck = (selector: DG.ColumnComboBox, col: string) =>
      selector.onChanged.pipe(filter((v) => v === col));

    const colSelection = 'There are a few ways to choose a column in a scatter plot. ' +
      'The easiest way to do this is to click on the column selector on the viewer. ' +
      'Alternatively, you can drag the column right from the spreadsheet, or from ' +
      'the column list (Windows | Columns, or Alt+C). Also, you can make this choice ' +
      'from the context panel on the right (Windows | Properties, or F4). ' +
      'Please try different ways in the next steps.';
    await this.action('Set X to HEIGHT', columnCheck(info.xColSelector, 'HEIGHT'), info.xColSelector.root, colSelection);
    await this.action('Set Y to WEIGHT',  columnCheck(info.yColSelector, 'WEIGHT'), info.yColSelector.root);
    await this.action('Set Size to AGE', columnCheck(info.sizeColSelector, 'AGE'), info.sizeColSelector.root);
    await this.action('Set Color to SEX', columnCheck(info.colorColSelector, 'SEX'), info.colorColSelector.root);

    const zoomDescription = 'To zoom in, hold the <b>Alt</b> key and drag a rectangle that you want to zoom in to.';
    await this.action('Zoom in', plot.onZoomed, null, zoomDescription);

    const zoomReset = 'As in most viewers, double-clicking on an empty area resets the view. ' +
      'Alternatively, this option is always available in the context menu.';
    await this.action('Double-click to unzoom', plot.onResetView, null, zoomReset);

    const currentRecord = 'Click on a point. Note that it becomes the current point in a spreadsheet, too.';
    await this.action('Click on a point', this.t!.onCurrentRowChanged, null, currentRecord);

    const selection = 'Select points by dragging a rectangle on a viewer while holding <b>Shift</b>. ' +
      'Note that the row selection is being reflected on most viewers, such as the spreadsheet. ' +
      'A number of points under the selection rectangle is shown right there.';
    await this.action('Select points', this.t!.onSelectionChanged, null, selection);

    const deselection = `Deselect some points by dragging a rectangle on a viewer while holding <b>Ctrl+Shift</b>.`;
    await this.action('Deselect points', this.t!.onSelectionChanged, null, deselection);
  }
}
