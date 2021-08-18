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

  protected async _run() {
    this.header.textContent = 'Scatter Plot';
    const plot = <DG.ScatterPlotViewer>(await this.openPlot(
      'scatter plot',
      (x) => x.type === DG.VIEWER.SCATTER_PLOT));
    const info = <{[key: string]: any}>plot.getInfo();
    const columnCheck = (selector: DG.ColumnComboBox, col: string) => selector.onChanged.pipe(filter((v) => v === col));

    this.describe("Now, let's customize it by choosing the columns of interest on corresponding axes. " +
      "The easiest way to do it is to click on the column selector on the viewer. Alternatively, " +
      "you can drag the column right from the spreadsheet, or from the column list " + 
      "(Windows | Columns, or Alt+C). Please try different ways to do it.");

    await this.action('Set X to HEIGHT', columnCheck(info.xColSelector, 'HEIGHT'), info.xColSelector.root);
    await this.action('Set Y to WEIGHT', columnCheck(info.yColSelector, 'WEIGHT'), info.yColSelector.root);
    await this.action('Set Size to AGE', columnCheck(info.sizeColSelector, 'AGE'), info.sizeColSelector.root);
    await this.action('Set Color to SEX', columnCheck(info.colorColSelector, 'SEX'), info.colorColSelector.root);
  }
}