import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import { fromEvent, interval, merge } from 'rxjs';
import $ from 'cash-dom';
import { waitForElementClick } from './utils';
import { getPlatform, Platform, platformKeyMap } from '../../shortcuts';

export class ViewersTutorial extends Tutorial {
  get name() { return 'Viewers'; }
  get description() {
    return 'Learn how to use different viewers together';
  }
  get steps() { return 20; }

  get icon() {
    return '📊👁️';
  }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers';
  platform: Platform = getPlatform();

  protected async _run() {
    this.showToolbox();
    this.header.textContent = this.name;
    this.title('Opening viewers');

    this.describe(`There are a few ways to add a viewer:
    <ol>
      <li>Toolbox (a panel on the left side)</li>
      <li>The top menu icon: <i class="grok-icon svg-icon svg-add-viewer"></i></li>
    </ol>
    `);
    this.describe('The icon opens a list of custom viewers, while ' +
    'the <b>Viewers</b> tab contains a standard set of visualizations.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    this.title('Selection and current records');

    const ribbonPanels = grok.shell.v.getRibbonPanels();
    const addViewerIcon = ribbonPanels[1][1];
    await this.action(
      'Click the Add viewer icon to open the gallery',
      waitForElementClick(addViewerIcon), addViewerIcon);

    const chartsSpanDiv = Array.from(document.querySelectorAll('.vg-tags .ui-div div')).find(div => {
      const span = div.querySelector('span');
      return span && span.textContent === 'Charts';
    }) as HTMLElement;
    await this.action('Click the "Charts" tag to filter the viewers', waitForElementClick(chartsSpanDiv), chartsSpanDiv);

    const radarViewerElement = this.getViewerCard('Radar');
    await this.action('Select the Radar viewer',
    waitForElementClick(radarViewerElement), radarViewerElement);

    const sunburstInfo = 'This time, find a viewer by searching instead of using tags. ' +
      '<b>Sunburst</b> shows hierarchical data as nested rings.';
    await this.action('Open the viewer gallery again', waitForElementClick(addViewerIcon), addViewerIcon, sunburstInfo);

    const searchInput = document.querySelector('.vg-controls input.ui-input-editor') as HTMLInputElement;
    await this.action('Type "Sunburst" in the search box',
      fromEvent(searchInput, 'input').pipe(filter(() => searchInput.value.toLowerCase().includes('sunburst'))),
      searchInput);

    const sunburstViewerElement = this.getViewerCard('Sunburst');
    await this.action('Select the Sunburst viewer',
      grok.events.onViewerAdded.pipe(filter((d: DG.EventData) => d.args.viewer.type === 'Sunburst')),
      sunburstViewerElement);

    const sp = await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    const hist = await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);
    const pie = await this.openPlot('pie chart', (x) => x.type === DG.VIEWER.PIE_CHART);

    const hover = 'Move the mouse over the histogram bins to see how the points ' +
    'that fall into that bin are reflected in other viewers. Similarly, hover the ' +
    'mouse over the pie chart segments or scatter plot data points.';
    await this.action('Hover over the histogram bins or scatter plot points',
      merge(this.t!.onMouseOverRowGroupChanged, this.t!.onMouseOverRowChanged), null, hover);

    const selection = 'Select points by dragging a rectangle on a viewer while holding <b>Shift</b>.';
    await this.action('Select points on the scatter plot', this.t!.onSelectionChanged, null, selection);

    const selectionSync = 'Note that the selection is synchronized between ' +
      'all viewers. When you select one of the bins on the histogram by clicking on it, ' +
      'you will see the corresponding records being highlighted on both scatter plot ' +
      'and grid. The same concept applies to the rest of the viewers, such as a pie chart ' +
      'or histogram. To select multiple data points, click on a segment while holding <b>Shift</b>. ' +
      `To deselect, hold <b>${platformKeyMap['Ctrl'][this.platform]}+Shift</b> while clicking. To invert, hold <b>${platformKeyMap['Ctrl'][this.platform]}</b> while clicking.`;
    await this.action('Select one of the bins on the histogram', this.t!.onSelectionChanged, null, selectionSync);

    const sunburstSelect = 'Click a <b>Sunburst</b> segment: every row under that branch is selected and synced to the other viewers.';
    await this.action('Click a Sunburst segment to select its rows', this.t!.onSelectionChanged, null, sunburstSelect);

    const currentRecord = 'Move the mouse over records on the scatter plot and grid, ' +
      'and note that the corresponding records are being highlighted in other viewers. ' +
      'Click on a point to make it current, and see how other viewers indicate where the current record is.';
    await this.action('Click on a point to set the current record', this.t!.onCurrentRowChanged, null, currentRecord);

    this.title('Properties');

    const spProperties = 'Make the scatter plot a current viewer by clicking on it, and then open its properties ' +
      '(press <b>F4</b> to bring out the context panel or click on the settings icon in the viewer header). There ' +
      'you can edit all properties of the viewer. Data-related properties are usually assembled on top under the ' +
      '<b>Data</b> category, while visual properties fall under various groups, such as <b>Colors</b>, ' +
      '<b>Markers</b> or <b>Axes</b>.';
    await this.action('Open the scatter plot\'s properties', interval(1000)
      .pipe(map((_) => grok.shell.o), filter((o) => o instanceof DG.ScatterPlotViewer)),
      null, spProperties);


    const initialProps = JSON.stringify((sp.getOptions() as {[key: string]: any }).look);
    await this.action('Change a few visual properties, e.g., the background color or marker size',
      interval(1000).pipe(map((_) => grok.shell.o), filter((o) => o instanceof DG.ScatterPlotViewer &&
      JSON.stringify((o.getOptions() as {[key: string]: any }).look) !== initialProps)));

    const cloneViewerInfo = 'Change some visual properties of the viewer and right-click on the scatter plot. ' +
      'In the context menu, select <b>General | Clone</b>. Note that the new viewer inherited all properties of the ' +
      'original viewer.<br> Close the new viewer by clicking on <b>"x"</b> in the top right corner of the header.';
    await this.contextMenuAction('Clone the scatter plot', 'Clone', null, cloneViewerInfo);

    this.title('Style');

    const stylePickInfo = 'You can apply the style of one viewer to another. ' +
      'To do that, right-click on the viewer and select <b>Pick Up / Apply | Pick Up</b>.';
    await this.contextMenuAction('Pick up the scatter plot\'s style', 'Pick Up', null, stylePickInfo);

    await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);

    const styleApplyInfo = 'To apply the style, choose <b>Pick Up / Apply | Apply</b> in the context menu ' +
      'of the new viewer. Depending on the situation, you might want to apply only visual or only ' +
      'data-related attributes; in this case, use <b>Apply Style Settings</b> or ' +
      '<b>Apply Data Settings</b>. Note that style settings can be applied even to ' +
      'viewers that have different source of data.';
    await this.contextMenuAction('Apply the style to the new viewer', 'Apply', null, styleApplyInfo);
  }

  private getViewerCard(name: string): HTMLElement | null {
    return Array.from(document.querySelectorAll<HTMLElement>('.viewer-gallery-root .d4-item-card.viewer-gallery'))
      .find(card => card.offsetParent !== null && card.querySelector('.card-label')?.textContent === name) ?? null;
  }
}
