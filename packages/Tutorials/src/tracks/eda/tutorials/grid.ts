import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';


export class GridTutorial extends Tutorial {
  get name(): string {
    return 'Grid Customization';
  }

  get description(): string {
    return '';
  }

  get steps(): number {
    return 10;
  }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers/grid';

  protected async _run(): Promise<void> {
    this.header.textContent = this.name;
    const view = grok.shell.tableView(this.t!.name);
    const grid = view.grid;
    this.title('Navigation');
    const navigationKeys = 'The current row in a grid is highlighted in green. To make a row current, ' +
      'click on it, or navigate up and down the grid using the arrow keys (<b>Up</b>, <b>Down</b>, ' +
      '<b>Left</b>, <b>Right</b>). To scroll faster, use <b>Page Up</b> and <b>Page Down</b>.';
    await this.action('Go a few rows down in the grid', this.t!.onCurrentRowChanged.pipe(filter(() =>
      this.t!.currentRow.idx > 0)), null, navigationKeys);

    await this.action('Jump back to the first row', this.t!.onCurrentRowChanged.pipe(filter(() =>
      this.t!.currentRow.idx === 0)), null, 'Press <b>Ctrl+Home</b> to return to the first row.');

    await this.action('Jump to the last row', this.t!.onCurrentRowChanged.pipe(filter(() =>
      this.t!.currentRow.idx === (this.t!.rowCount - 1))), null, 'Press <b>Ctrl+End</b> to go to the last row.');

    const columnNames = this.t!.columns.names();
    const columnNumber = this.t!.columns.length;
    await this.action('Jump to the last column', this.t!.onCurrentColChanged.pipe(filter(() =>
      this.t!.currentCol.name === columnNames[columnNumber - 1])), null, 'Press <b>End</b> to go to the last column.');
    await this.action('Jump to the first column', this.t!.onCurrentColChanged.pipe(filter(() =>
    this.t!.currentCol.name === columnNames[0])), null, 'Press <b>Home</b> to return to the first column.');

    // Edit values
    await this.action('Click the <b>+</b> icon', this.t!.onRowsAdded, $('.d4-grid-add-row-icon')[0], '');
    this.title('Selection');
    // Selection of rows/columns
    // Add/remove rows/columns
    this.title('Sorting');
    // Sort
    // Resize/reorder
    this.title('Formatting');
    // Formatting of dates and numbers
    this.title('Color coding');
    // Color coding
    this.title('Summary columns');
    // Summary columns
    // Hide columns
  }
}
