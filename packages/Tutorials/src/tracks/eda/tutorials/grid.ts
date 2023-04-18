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

    this.title('Editing');

    await this.action('Add a new row by clicking the "+" icon at the last row', this.t!.onRowsAdded,
      $('.d4-grid-add-row-icon')[0], 'There are several ways to add a row to a spreadsheet: clicking the <b>+</b> ' +
      'icon (requires "Allow Edit" property set to true), editing the last row (if "Add New Row On Last Row Edit" ' +
      'is enabled), selecting <b>Edit > Add Rows...</b> in the top menu (this way you can add a specific number of ' +
      'rows at a specific position in the grid).');

    const lastRowIdx = this.t!.rowCount - 1;
    await this.action(`Enter "X0273T51080200024" as USUBJID in the new row (#${lastRowIdx + 1})`, this.t!.onValuesChanged.pipe(filter(() =>
      this.t!.cell(lastRowIdx, 'USUBJID').value === 'X0273T51080200024')), null, 'Once you finish editing, ' +
      'one more row will appear (due to the above mentioned "Add New Row On Last Row Edit" grid property).');

    await this.action(`Set AGE to "37" in the row you edited (#${lastRowIdx + 1})`, this.t!.onValuesChanged.pipe(
      filter(() => this.t!.cell(lastRowIdx, 'AGE').value === 37)));

    await this.action(`Set SEX to "F" in the row you edited (#${lastRowIdx + 1})`, this.t!.onValuesChanged.pipe(
      filter(() => this.t!.cell(lastRowIdx, 'SEX').value === 'F')));

    await this.action(`Copy the value from row #${lastRowIdx} to row #${lastRowIdx + 1} for the RACE column`,
      this.t!.onValuesChanged.pipe(filter(() => this.t!.cell(lastRowIdx, 'RACE').value ===
      this.t!.cell(lastRowIdx - 1, 'RACE').value)), null, 'Use <b>Ctrl+C</b> to copy the value of the ' +
      'current cell and <b>Ctrl+V</b> to paste into the current cell.');

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
