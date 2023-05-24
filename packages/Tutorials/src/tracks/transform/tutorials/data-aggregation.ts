import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import $ from 'cash-dom';
import { interval } from 'rxjs';
import { filter } from 'rxjs/operators';



export class AggregationTutorial extends Tutorial {
  get name(): string {
    return 'Data Aggregation';
  }
  get description(): string {
    return 'Learn different ways of data aggregation and pivoting';
  }
  get steps(): number {
    return 10;
  }

  helpUrl: string = 'https://datagrok.ai/help/transform/aggregate-rows';

  protected async _run() {
    grok.shell.windows.showToolbox = false;
    this.header.textContent = this.name;
    this.title('Aggregation');

    this.describe('<b>Aggregation Editor</b> is a visual tool for summarizing and pivoting ' +
      'table data. It comes with a wide set of aggregation functions. We will learn how to ' +
      'aggregate data on a demographics dataset.');

    await this.action('Open Aggregation Editor', grok.functions.onAfterRunAction.pipe(
      filter((call) => call.func.name === 'CmdAggregateRows')), this.getMenuItem('Data'),
      'Select <b>Data > Aggregate Rows</b> in the top menu, or press <b>Alt+A</b>.');

    this.describe('The aggregation editor consists of several components: the section on ' +
      'top contains aggregation parameters; the spreadsheet at the bottom shows the result, ' +
      'which is calculated interactively as the input changes; the column manager on the left ' +
      'lists all columns of the source dataframe. You can see column stats in the tooltip and ' +
      'drag relevant columns to the editor on the right in order to aggregate or pivot by them.');

    const colTagSelector = 'div.d4-tag-editor.d4-pivot-column-tags span.d4-tag';
    const findColTag = (root: HTMLElement, colName: string, condition: (() => boolean) | null = null) =>
      interval(1000).pipe(filter(() => {
        const tag = $(root).find(colTagSelector)
          .filter((idx, el) => el.textContent?.toUpperCase() === colName.toUpperCase())[0];
        return tag != null && (condition === null ? true : condition());
    }));

    const groupByRoot = $('.grok-pivot-column-panel').filter((_, el) =>
      el.textContent?.startsWith('Group by') === true)[0]!;

    const groupByCol1 = 'RACE';
    await this.action(`Group rows by column "${groupByCol1}"`, findColTag(groupByRoot, groupByCol1),
      groupByRoot.querySelector('.grok-pivot-column-tags-plus') as HTMLElement, 'To group rows, put ' +
      'the corresponding column in the "Group by" field. Click on the <b>"+"</b> sign next to the ' +
      'field header and select the column.');

    const groupByCol2 = 'SEX';
    await this.action(`Group rows by column "${groupByCol2}"`, findColTag(groupByRoot, groupByCol2),
      null, 'Another way to add a column to group by is to drag it from the column list on the right. ' +
      'Add the second column to the grouping list.');

    this.title('Pivoting');

    this.describe('Pivoting is a way to transform data to group values belonging to the same ' +
      'category in columns instead of rows. Let\'s pivot our dataset by the patient\'s condition.');

    const pivotRoot = $('.grok-pivot-column-panel').filter((_, el) =>
      el.textContent?.startsWith('Pivot') === true)[0]!;

    const pivotCol = 'DIS_POP';
    await this.action(`Pivot data by column "${pivotCol}"`, findColTag(pivotRoot, pivotCol),
      pivotRoot.querySelector('.grok-pivot-column-tags-plus') as HTMLElement);

    this.title('Aggregation');

    const aggRoot = $('.grok-pivot-column-panel').filter((_, el) =>
      el.textContent?.startsWith('Aggregate') === true)[0]!;

    await this.action('Leave only the "avg(AGE)" aggregation',
      findColTag(aggRoot, 'avg(AGE)', () => $(aggRoot).find(colTagSelector).length === 1), null,
      'Initially, the editor shows the average values for the "AGE" and "HEIGHT" columns. To keep ' +
      'only one aggregation, right-click this aggregation and select <b>Remove others</b> in the ' +
      'context menu.');

    await this.action('Change a column to "WEIGHT"', findColTag(aggRoot, 'avg(WEIGHT)', () =>
      $(aggRoot).find(colTagSelector).length === 1), null, 'You can change the aggregation column ' +
      'and function from the context menu. For example, if you are adding multiple columns using the ' +
      'same aggregation function, you can set it as default by pressing the "+" sign and choosing ' +
      'it under the "Aggregation" submenu. Let\'s calculate the average weight for each patient ' +
      'group instead of age. Right-click the aggregation field and select <b>Column > WEIGHT</b>.');

    this.title('Interactivity');

    this.describe('To make data exploration easier, the aggregated table is linked with the source ' +
      'table. Whenever a row is selected in the aggregated table, corresponding rows get selected ' +
      'in the source table. When current record changes, the source table gets filtered to show ' +
      'only rows associated with the current record.');

    await this.action('Select rows in the source table with values of the first aggregated row',
      this.t!.selection.onChanged.pipe(filter(() => this.t!.selection.trueCount === 37)), null,
      'Click on the first row in the aggregated table while holding <b>Shift</b>. This way you ' +
      'will select all the corresponding rows in the source table (the values are "Asian, F").');

    await this.action('Click on the last row in the aggregated table to filter by it',
      this.t!.filter.onChanged.pipe(filter(() => this.t!.filter.trueCount === 75)), null,
      'The filter should be based on the last row in the aggregated table (the values are ' +
      '"Other, M").');

    this.title('History');

    const serializedParams = '[{"#type":"GroupAggregation","aggType":"key","colName":"RACE"},' +
      '{"#type":"GroupAggregation","aggType":"key","colName":"SEX"},{"#type":"GroupAggregation",' +
      '"aggType":"pivot","colName":"DIS_POP"},{"#type":"GroupAggregation","aggType":"avg","colName":"WEIGHT"}]';

    const historyStr = window.localStorage['grok-aggregation-history'];
    const initialParamsLen = historyStr ? JSON.parse(historyStr).length : 0;

    await this.action('Save parameters', interval(1000).pipe(filter(() => {
        const historyStr = window.localStorage['grok-aggregation-history'];
        if (!historyStr)
          return false;
        const history = JSON.parse(historyStr);
        return history.length > initialParamsLen &&
          JSON.stringify(history[history.length - 1]) === serializedParams;
      })), $('i.grok-icon.fa-history.d4-command-bar-icon')[0],
      'Click on the history icon and select <b>Save parameters</b> from the menu. Note that parameters are ' +
      'also saved automatically when you click "OK" to add the aggregated dataframe to the workspace. This ' +
      'can be useful if you choose to reset the entered parameters. To return previously used parameters, ' +
      'simply select them when clicking the history icon.');
  }
}
