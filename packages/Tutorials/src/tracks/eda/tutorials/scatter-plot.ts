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
  get steps(){ return 10} //set thew number of steps
    
  protected async _run() {
    this.header.textContent = 'Scatter Plot';

    this.describe('A scatter plot (also called a scatter graph, scatter chart, scattergram, or scatter diagram) is a type of plot or mathematical diagram using Cartesian coordinates to display values for typically two variables for a set of data. If the points are color-coded you can increase the number of displayed variables to three. The data is displayed as a collection of points, each having the value of one variable determining the position on the horizontal axis and the value of the other variable determining the position on the vertical axis.');

    const plot = <DG.ScatterPlotViewer>(await this.openPlot(
      'scatter plot','',
      (x) => x.type === DG.VIEWER.SCATTER_PLOT));
    const info = <{ [key: string]: any }>plot.getInfo();
    const columnCheck = (selector: DG.ColumnComboBox, col: string) => selector.onChanged.pipe(filter((v) => v === col));

    await this.action('Set X to HEIGHT', 
    "Now, let's customize it by choosing the columns of interest on corresponding axes. The easiest way to do it is to click on the column selector on the viewer. Alternatively, you can drag the column right from the spreadsheet, or from the column list (Windows | Columns, or Alt+C). Please try different ways to do it.", columnCheck(info.xColSelector, 'HEIGHT'), info.xColSelector.root);
    await this.action('Set Y to WEIGHT', '',  columnCheck(info.yColSelector, 'WEIGHT'), info.yColSelector.root);
    await this.action('Set Size to AGE', '', columnCheck(info.sizeColSelector, 'AGE'), info.sizeColSelector.root);
    await this.action('Set Color to SEX', '', columnCheck(info.colorColSelector, 'SEX'), info.colorColSelector.root);

    await this.action('Zoom in','To zoom in, hold Alt key, and drag a rectangle that you want to zoom in to. Try it now.', plot.onZoomed);

    await this.action('Double-click to unzoom', 'As in most viewers, double-clicking on an empty area resets the view. Alternatively, this option is always available in the context menu.', plot.onResetView);

    await this.action('Click on a point', 'Click on a point. Note that it becomes the current point in a spreadsheet, too.',this.t.onCurrentRowChanged);

    await this.action('Select points', 'Select points by dragging a rectangle on a viewer while holding Shift. Note that the row selection is being reflected on most viewers, such as the spreadsheet. A number of points under the selection rectangle is shown right there.',this.t.onSelectionChanged);

    await this.action('Deselect points',`Deselect some points by dragging a rectangle on a viewer while holding <b>Ctrl+Shift.</b>`, this.t.onSelectionChanged);
  }
}
