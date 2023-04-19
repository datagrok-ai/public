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
    await this.action(`Enter "X0273T51080200024" as USUBJID in the new row (#${lastRowIdx + 1})`,
      this.t!.onValuesChanged.pipe(filter(() => this.t!.cell(lastRowIdx, 'USUBJID').value === 'X0273T51080200024')),
      null, 'Once you finish editing, one more row will appear (due to the above mentioned "Add New Row On Last ' +
      'Row Edit" grid property).');

    await this.action(`Set AGE to "37" in the row you edited (#${lastRowIdx + 1})`, this.t!.onValuesChanged.pipe(
      filter(() => this.t!.cell(lastRowIdx, 'AGE').value === 37)));

    await this.action(`Set SEX to "F" in the row you edited (#${lastRowIdx + 1})`, this.t!.onValuesChanged.pipe(
      filter(() => this.t!.cell(lastRowIdx, 'SEX').value === 'F')));

    await this.action(`Copy the value from row #${lastRowIdx} to row #${lastRowIdx + 1} for the RACE column`,
      this.t!.onValuesChanged.pipe(filter(() => this.t!.cell(lastRowIdx, 'RACE').value ===
      this.t!.cell(lastRowIdx - 1, 'RACE').value)), null, 'Use <b>Ctrl+C</b> to copy the value of the ' +
      'current cell and <b>Ctrl+V</b> to paste the copied value into the current cell.');

    this.title('Selection');

    await this.action(`Select the edited row by clicking on its number (#${lastRowIdx + 1})`,
      this.t!.onSelectionChanged.pipe(filter(() => this.t!.selection.trueCount === 1 && this.t!.selection.get(lastRowIdx))),
      null, 'Grid shares row selection with other viewers, which will be illustrated in the upcoming tutorials. ' +
      'Selected data is highlighted in orange. You can extend selection by holding <b>Shift</b> and dragging the mouse.');

    await this.action('Remove selection by pressing "Esc"', this.t!.onSelectionChanged.pipe(filter(() =>
      !this.t!.selection.anyTrue)), null, 'There are multiple ways to deselect all or certain rows, hitting ' +
      '<b>Esc</b> being the simplest. It removes selection entirely, both from rows and columns.');

    await this.action('Select all rows with "Ctrl+A"', this.t!.onSelectionChanged.pipe(filter(() =>
      !this.t!.selection.anyFalse)), null, 'There is a key combination to revert this action. You can use ' +
      '<b>Ctrl+Shift+A</b> to deselect all rows. It differs from <b>Esc</b>, since it affects only row selection.');

    const rowCount = this.t!.rowCount;
    const selectedIndexes = [0, 1, 2, 3, 4, rowCount - 5, rowCount - 4, rowCount - 3, rowCount - 2, rowCount - 1];
    await this.action(`Select the first five and the last five table rows (#1 - #5 and #${rowCount - 4} - #${rowCount})`,
      this.t!.onSelectionChanged.pipe(filter(() => this.t!.selection.trueCount === 10 && selectedIndexes.every((idx, i) =>
      this.t!.selection.getSelectedIndexes()[i] === idx))), null, 'Make use of the previously learned shortcuts. ' +
      '<b>Shift+Drag</b> selects rows, <b>Ctrl+Shift+Drag</b> deselects rows, <b>Ctrl+Click</b> toggles selected state. ' +
      '<b>Ctrl+Shift+Home</b> selects rows above the current row, while <b>Ctrl+Shift+End</b> selects rows below the ' +
      'current one (inclusively).');

    await this.action('Clear the selection', this.t!.onSelectionChanged.pipe(filter(() => !this.t!.selection.anyTrue)),
      null, 'Hit either <b>Esc</b> or <b>Ctrl+Shift+A</b>.');

    grok.shell.windows.showContextPanel = true;
    const noneSeverityRowCount = 3302;
    const severityColumn = this.t!.getCol('SEVERITY');
    await this.action('Find the "SEVERITY" column and select all rows with the "None" value', this.t!.onSelectionChanged.pipe(
      filter(() => this.t!.selection.trueCount === noneSeverityRowCount)), null, 'First, find any cell with "None" severity ' +
      'and click to make it current. Then press <b>Ctrl+Enter</b>, it will select all rows containing the same value.');

    await this.action(`Delete these rows (${noneSeverityRowCount})`, this.t!.onRowsRemoved.pipe(filter(() =>
      this.t!.rowCount === (rowCount - noneSeverityRowCount) && !severityColumn.categories.includes('None'))),
      null, 'In the context panel on the right, you can find the list of actions for selected rows. One of ' +
      'them is "Delete Rows". Pick this option to refine your dataset.');
    // Selection of columns
    // Add and remove columns
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
