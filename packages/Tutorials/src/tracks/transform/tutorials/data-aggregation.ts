import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import $ from 'cash-dom';
import { interval } from 'rxjs';
import { filter } from 'rxjs/operators';



export class AggregationTutorial extends Tutorial {
  get name() {
    return 'Data Aggregation';
  }
  get description() {
    return 'Learn different ways of data aggregation and pivoting';
  }
  get steps() {
    return 10;
  }

  helpUrl: string = 'https://datagrok.ai/help/transform/aggregate-rows';

  protected async _run() {
    this.title('Aggregation');

    await this.action('Open Aggregation Editor', grok.functions.onAfterRunAction.pipe(
      filter((call) => call.func.name === 'CmdAggregateRows')), this.getMenuItem('Data'),
      'Select <b>Data > Aggregate Rows</b> in the top menu, or press <b>Alt+A</b>.');

    this.describe('The aggregation editor consists of several components: the section on ' +
      'top contains aggregation parameters; the spreadsheet at the bottom shows the result, ' +
      'which is calculated interactively as the input changes; the column manager on the left ' +
      'lists all columns of the source dataframe. You can see column stats in the tooltip and ' +
      'drag relevant columns to the editor on the right in order to aggregate or pivot by them.');

    const groupByCol1 = 'RACE';
    await this.action(`Group rows by "${groupByCol1}"`, interval(1000).pipe(filter(() =>
      $('div.d4-tag-editor.d4-pivot-column-tags span.d4-tag').filter((idx, el) =>
      el.textContent?.toUpperCase() === groupByCol1)[0] != null)),
      $('.grok-pivot-column-panel').filter((idx, el) => el.textContent?.startsWith('Group by') ===
      true).find('.grok-pivot-column-tags-plus')[0], '');
  }
}
