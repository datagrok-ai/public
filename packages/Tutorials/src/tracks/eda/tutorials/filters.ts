import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';
import wu from 'wu';


export class FiltersTutorial extends Tutorial {
  get name() {
    return 'Filters';
  }
  get description() {
    return 'A set of controls for quick filtering, selection, ' +
      'and visual assessment of column values';
  }

  get steps() { return 18; }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers/filters';
  
  protected async _run() {
    this.header.textContent = this.name;
    this.describe('Dynamic filtering is an important concept in exploratory data analysis, and ' +
      'our platform makes it as powerful and easy to use as possible. Let\'s start with opening ' +
      'a couple of regular viewers, so that effects of filtering would be immediately visible.');

    this.describe(ui.link('More about ' + this.name, this.helpUrl).outerHTML);

    await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);

    const filterDescription = 'While appearing very simple, this is a powerful tool. ' +
      'It combines multiple indicators along with a few controls that let you easily ' +
      `modify current filters and selection. Let's start with filtering by a specific category.`;
    const filters = await this.openPlot('filters', (x) => x.type === DG.VIEWER.FILTERS, filterDescription);

    this.title('Categorical filters');

    const catFilterInfo = 'Only rows that have this value in the <b>DIS_POP</b> column ' +
      'will pass the filter. You can see the number of filtered rows and the applied filter ' +
      'description in the property panel on the right. In categorical filters, the blue ' +
      'background, as well as the counts on the right, indicates the ratio of rows for that ' +
      'category that pass the current filter. We have just filtered the dataset by disease ' +
      'type, so naturally the rest of categories have zero row counts. Note that this is ' +
      'not the case for other columns, for example, based on the state of the <b>SEX</b> ' +
      'filter, we can immediately tell that within patients with ankylosing spondylitis ' +
      '(AS), there is a higher proportion of men.';
    await this.action('Click on the "AS" label within the "DIS_POP" filter',
      this.t!.onFilterChanged.pipe(filter(() => {
        const filters = this.t!.rows.filters;
        return filters.length === 1 && filters.get(0) === 'DIS_POP: AS';
      })), null, catFilterInfo);

    const keysToSwitch = 'A convenient way to quickly browse data by categories is to click ' +
      'on any of them, and then use the up and down cursor keys to move between the categories ' +
      'in the filter. Pay attention to the row counts, it seems that rheumatoid arthritis (RA) ' +
      'occurs quite often.';
    await this.action('Filter the dataset by the most common disease',
    this.t!.onFilterChanged.pipe(filter(() => {
      const filters = this.t!.rows.filters;
      return filters.length === 1 && filters.get(0) === 'DIS_POP: RA';
    })), null, keysToSwitch);

    const indicatorInfo = 'To filter by multiple categories, either check each label manually ' +
      'or right-click the indicator beside the column name and choose <b>Invert all</b>.';
    await this.action('Invert the category set for the "DIS_POP" filter',
      this.t!.onFilterChanged.pipe(filter(() => {
        const filters = this.t!.rows.filters;
        return filters.length === 1 &&
          filters.get(0) === 'DIS_POP: AS, Indigestion, PsA, Psoriasis, UC';
      })), $('div.d4-flex-row.d4-filter-header').filter((idx, el) => $(el)
        .find('label.d4-filter-column-name')[0]?.textContent === 'DIS_POP')
        .find('div.d4-filter-indicator')[0],
      indicatorInfo);

    const rowCountSelect = 'When you click on a row count, the corresponding rows ' +
      'get selected, taking into account the current filter. They are highlighted ' +
      'in orange both in filters and other viewers.';
    await this.action('Click on a non-empty row count', this.t!.onSelectionChanged, null, rowCountSelect);

    await this.action('Filter the dataset to only females of Asian or Black origin',
      this.t!.onFilterChanged.pipe(filter(() => {
        let rMatch, sMatch;
        rMatch = sMatch = false;
        wu(this.t!.rows.filters).forEach((f) => {
          if (f === 'RACE: Asian, Black') {
            rMatch = true;
          } else if (f === 'SEX: F') {
            sMatch = true;
          }
        });
        return rMatch && sMatch;
      })),
      $('label.d4-filter-column-name').filter((idx, el) =>
        el.textContent === 'SEX' || el.textContent === 'RACE').get(),
      'Combine two filters: <b>SEX</b> and <b>RACE</b>.');

    await this.action('Reset the filter',
      this.t!.onFilterChanged.pipe(filter((_) => this.t!.filter.trueCount === this.t!.rowCount)),
      $(filters.root).find('i.grok-icon.fa-sync')[0],
      'Press <b>Esc</b> or click on <i class="grok-icon fal fa-sync"></i> at the top of the filter panel.');

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
      this.t!.onFilterChanged.pipe(filter(() => wu(this.t!.rows.filters).some((s) => s === 'AGE: [40,60]'))),
      $('div.d4-flex-row.d4-filter-header').filter((idx, el) => $(el)
        .find('label.d4-filter-column-name')[0]?.textContent === 'AGE')
        .find('i.grok-icon.fa-keyboard')[0],
      rangeInputInfo);

    this.title('Saving filter state');

    await this.contextMenuAction('Save the current filter as a column', 'Filter to Column...', null,
      'To save the current filter state into a boolean column, use the <b>Filter to column...</b>' +
      ' command in the context menu. A new column will be added to the dataframe. By default, ' +
      'its name corresponds to the applied filters, for example, <b>AGE: [40,60]</b>, but you ' +
      'can change it from the dialog when saving.');

    await this.contextMenuAction('Save the filter configuration as "AGE: [40,60]"', 'Save...', null,
      'Choose <b>Save or apply | Save...</b> in the context menu ' +
      'and keep the default name <b>AGE: [40,60]</b>.');

    await this.action('Reset the filter',
      this.t!.onFilterChanged.pipe(filter((_) => this.t!.filter.trueCount === this.t!.rowCount)),
      $(filters.root).find('i.grok-icon.fa-sync')[0],
      'Press <b>Esc</b> or click on <i class="grok-icon fal fa-sync"></i> at the top of the filter panel.');

    await this.contextMenuAction('Restore the filter state', 'AGE: [40,60]', null,
      'Find the saved filter state in <b>Save or apply</b> and click on its name. You should ' +
      'see that the dataset is filtered to patients between the ages of 40 and 60 again.');
  }
}
