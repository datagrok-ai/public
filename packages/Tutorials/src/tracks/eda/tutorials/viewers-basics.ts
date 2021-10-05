import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';
import { interval } from 'rxjs';


export class ViewersTutorial extends Tutorial {
  get name() { return 'Viewers'; }
  get description() {
    return 'Learn how to use different viewers together';
  }
  get steps() { return 9; }
  
  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers';

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

    this.describe(String(ui.link('More about '+this.name, this.helpUrl).outerHTML));

    this.title('Selection and the current record');

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
    await this.action('Click on a point to set the current record', this.t!.onCurrentRowChanged, null, currentRecord);

    this.title('Properties');

    const spProperties = 'Make the scatter plot a current viewer by clicking on it, and then open its properties ' +
      '(press F4 to bring out the property panel or click on the settings icon in the viewer header). There you ' +
      'can edit all properties of the viewer. Data-related properties are usually assembled on top under the ' +
      '"Data" category, while visual properties fall under various groups, such as "Colors", "Markers" or "Axes".' +
      '<br> Go ahead and change some of the appearance properties of the scatter plot, such as the background color.'
    await this.action('Open the scatter plot\'s properties', interval(1000)
      .pipe(map((_) => grok.shell.o), filter((o) => o instanceof DG.Viewer)),
      null, spProperties);

    const cloneViewerInfo = 'Change some visual properties of the viewer and right-click on the scatter plot. ' +
      'In the context menu, select "General > Clone". Note that the new viewer inherited all properties of the ' +
      'original viewer. <br> Close the new viewer by clicking on "x" in the top right corner of the header.';
    await this.contextMenuAction('Clone the scatter plot', 'Clone', null, cloneViewerInfo);
  }
}
