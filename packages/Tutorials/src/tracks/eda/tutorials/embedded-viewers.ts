import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter, map } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import { interval } from 'rxjs';


export class EmbeddedViewersTutorial extends Tutorial {
  get name(): string {
    return 'Embedded Viewers';
  }
  get description(): string {
    return 'Find out how to use viewers in tooltips and inside other viewers';
  }
  get steps(): number {
    return 10;
  }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers';

  protected async _run(): Promise<void> {
    this.header.textContent = this.name;
    this.describe('Viewers can be placed in tooltips or be part of other viewers. ' +
      'In this tutorial, we will learn how to change the standard tooltip template ' +
      'and work with a Trellis plot that incorporates other viewers.');
    this.describe(ui.link('More about viewers', this.helpUrl).outerHTML);

    this.title('Viewer tooltips');

    const sp = await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    const hist = await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);
    const pie = await this.openPlot('pie chart', (x) => x.type === DG.VIEWER.PIE_CHART);

    await this.contextMenuAction('Set the scatter plot as a tooltip viewer',
      'Use as group tooltip', null, 'Each viewer can be used as a template for the mouse-over ' +
      'rows to be visualized in a tooltip. To understand how it works, make the scatter plot ' +
      'a tooltip viewer by selecting <b>Tooltip | Use as Group Tooltip</b> from the context menu.');

    await this.action('Hover over the histogram bins or pie chart segments',
      this.t!.onMouseOverRowGroupChanged, null,
      'The rows represented by each bin will be visualized in the tooltip as the scatter plot. ' +
      'Try doing the same with the pie chart segments to see scatter plots for different row groups.');

    await this.contextMenuAction('Reset the tooltip', 'Remove Group Tooltip', null,
      'Choose the option <b>Tooltip | Remove Group Tooltip</b> in the context menu.');

    this.title('Trellis plot');

    await this.contextMenuAction('Open a Trellis plot from the pie chart\'s context menu',
      'Use in Trellis', null, 'The most common viewers can be put on a trellis plot by choosing ' +
      '<b>General | Use in Trellis</b>. You can select columns by which the data should be split. ' +
      'Each cell will only show rows that belong to the corresponding categories.');

    let trellis: DG.Viewer;
    for (const v of (<DG.TableView>grok.shell.v).viewers) {
      if (v.type === DG.VIEWER.TRELLIS_PLOT) {
        trellis = v;
      }
    }

    await this.action('Set a scatter plot as an inner viewer', interval(1000).pipe(
        map((_) => trellis.props.viewerType),
        filter((t: string) => t === DG.VIEWER.SCATTER_PLOT)),
      $(trellis!.root).find('.d4-combo-popup')[0], 'This time, use the viewer type selector ' +
      'in the opened Trellis plot. You can also set it from the context panel.');

    await this.action('Open the inner plot properties and set Color to AGE',
      interval(1000).pipe(
        map((_) => (trellis.getOptions() as {[key: string]: any }).look.innerViewerLook.colorColumnName),
        filter((name: string | undefined) => name === 'AGE')),
      $(trellis!.root).find('.grok-font-icon-settings.d4-viewer-icon')[0],
      'Click on the gear icon to edit the scatter plot properties. ' +
      'Find the <b>Color</b> property within the <b>Color</b> section.');
  }
}
