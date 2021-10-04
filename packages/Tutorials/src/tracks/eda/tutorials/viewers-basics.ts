import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';
import { Observable } from 'rxjs';


export class ViewersTutorial extends Tutorial {
  get name() { return 'Viewers'; }
  get description() {
    return 'Learn how to use different viewers together';
  }
  get steps() { return 6; }
    
  protected async _run() {
    this.title('Opening viewers');

    this.describe(`There are a few ways to add a viewer:
    <ol>
      <li>Toolbox (a panel on the left side)</li>
      <li>The top menu icon: <i class="grok-icon svg-icon svg-add-viewer"></i></li>
    </ol>
    `);
    this.describe('The icon opens a list of custom viewers, while ' +
    'the "Viewers" tab contains a standard set of visualizations.');

    this.describe("Let's start by opening some viewers.");

    const sp = await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    const hist = await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);
    const pie = await this.openPlot('pie chart', (x) => x.type === DG.VIEWER.PIE_CHART);

    const selection = 'Select points by dragging a rectangle on a viewer while holding <b>Shift</b>.';
    await this.action('Select points on the scatter plot', this.t!.onSelectionChanged, null, selection);

    const selectionSync = 'Move the mouse over histogram bins to see how the points ' +
      'that fall into that bin are reflected in other viewers. Similarly, hover the ' +
      'mouse over pie chart segments. Note that the selection is synchronized between ' +
      'all viewers. When you select one of the bins on the histogram by clicking on it, ' +
      'you will see the corresponding records being highlighted on both scatter plot ' +
      'and grid. The same concept applies to the rest of the viewers, such as a pie chart ' +
      'or histogram. To select multiple data points, click on a segment while holding <b>Shift</b>. ' +
      'To deselect, hold <b>Ctrl+Shift</b> while clicking. To invert, hold <b>Ctrl</b> while clicking.';
    await this.action('Select one of the bins on the histogram', this.t!.onSelectionChanged, null, selectionSync);

    const currentRecord = 'Move the mouse over records on the scatter plot and grid, ' +
      'and note that the corresponding records are being highlighted in other viewers. ' +
      'Click on a point to make it current, and see how other viewers indicate where the current record is.';
    await this.action('Click on a point to set the current record.', this.t.onCurrentRowChanged, null, currentRecord);

    await this.action('Now, select one of the bins on the histogram', `Move the mouse over histogram bins to see how the points that fall into that bin are reflected in other viewers. Similarly, hover the mouse over pie chart segments.\n\nNote that the selection is synchronized between all viewers. Now, select one of the bins on the histogram by clicking on it, and see the corresponding records being highlighted on both scatter plot and grid.\n\nThe same concept applies to the rest of the viewers, such as a pie chart or histogram. To select multiple data points, click on a segment while holding <b>Shift</b>. To unselect, hold <b>Ctrl+Shift</b> while clicking. To invert, hold <b>Ctrl</b> while clicking.`, this.t.onSelectionChanged);
    
    await this.action('Click on a point to set the current record.', 'Move the mouse over records on the scatter plot and grid, and note that the corresponding records are being highlighted in other viewers. Click on a point to make it current, and see how other viewers indicate where the current record is.' , this.t.onCurrentRowChanged);

  }
}
