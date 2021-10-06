import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class FiltersTutorial extends Tutorial {
  get name() {
    return 'Filters';
  }
  get description() {
    return 'A set of controls for quick filtering, selection, ' +
      'and visual assessment of column values';
  }

  get steps() { return 8; }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers/filters';
  
  protected async _run() {
    this.header.textContent = 'Filters';
    this.describe('Dynamic filtering is an important concept in exploratory data analysis, and ' +
      'our platform makes it as powerful and easy to use as possible. Let\'s start with opening ' +
      'a couple of regular viewers, so that effects of filtering would be immediately visible.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);

    const filterDescription = 'While appearing very simple, this is a very powerful tool. ' +
      'It combines multiple indicators along with a few controls that let you easily ' +
      `modify current filters and selection. Let's start with filtering by a specific category.`;
    const filters = await this.openPlot('filters', (x) => x.type === DG.VIEWER.FILTERS, filterDescription);

    // this.title('Categorical filters');
    this.title('Numerical filters');

    const hoverInfo = 'Move the mouse over the histogram bins in a numeric filter ' +
      'or in the viewer. You will see how many rows fall into each bin\'s range.';
    await this.action('Hover over the histogram bins', this.t!.onMouseOverRowGroupChanged, null, hoverInfo);

    const selectionInfo = 'Click on a bin to select it. To select multiple bins, hold <b>Shift</b> and ' +
      'either pick bins one by one or drag a rectangle to form a group. If you hold <b>Ctrl</b> while ' +
      'clicking, you will toggle the bin\'s selection. Note that other filters reflect the proportion ' +
      'of the selected rows.';
    await this.action('Select one of the histogram bins', this.t!.onSelectionChanged, null, selectionInfo);

    const indicatorsInfo = 'Current and mouse-over records are shown below the histogram ' +
      'bins as green and gray circles. These indicators can be used for quick data profiling.';
    await this.action('Change the current row in the spreadsheet', this.t!.onCurrentRowChanged, null, indicatorsInfo);

    const rangeInputInfo = 'When the mouse is over a histogram, a range slider appears at the bottom. ' +
      'By dragging the handles at the edges of the slider or panning it, you can define the range of ' +
      'the values that pass filter. For more accurate results, use the <i class="grok-icon fal fa-keyboard"></i> ' +
      'icon in the numeric filter header. It toggles the range inputs.';
    await this.action('Find records for people aged 40 to 60',
      this.t!.onFilterChanged.pipe(filter(() => this.t!.rows.filters.some((s) => s === 'AGE: [40,60]'))),
      $('div.d4-flex-row.d4-filter-header').filter((idx, el) => $(el)
        .find('label.d4-filter-column-name')[0]?.textContent === 'AGE')
        .find('i.grok-icon.fa-keyboard')[0],
      rangeInputInfo); 
  }
}
