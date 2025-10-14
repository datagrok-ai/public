import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import { interval } from 'rxjs';
import { filter } from 'rxjs/operators';
import { Tutorial } from '@datagrok-libraries/tutorials/src/tutorial';
import { Platform, getPlatform, platformKeyMap } from '../../shortcuts';


export class GridTutorial extends Tutorial {
  get name(): string {
    return 'Grid Customization';
  }

  get description(): string {
    return 'The main viewer for interactive exploration of tables';
  }

  get steps(): number {
    return 31;
  }

  helpUrl: string = 'https://datagrok.ai/help/visualize/viewers/grid';
  platform: Platform = getPlatform();

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
      this.t!.currentRow.idx === 0)), null, `Press <b>${platformKeyMap['Ctrl + Home'][this.platform]}</b> to return to the first row.`);

    await this.action('Jump to the last row', this.t!.onCurrentRowChanged.pipe(filter(() =>
      this.t!.currentRow.idx === (this.t!.rowCount - 1))), null, `Press <b>${platformKeyMap['Ctrl + End'][this.platform]}</b> to go to the last row.`);

    const columnNames = this.t!.columns.names();
    const columnNumber = this.t!.columns.length;
    await this.action('Jump to the last column', this.t!.onCurrentColChanged.pipe(filter(() =>
      this.t!.currentCol.name === columnNames[columnNumber - 1])), null, `Press <b>${platformKeyMap['End'][this.platform]}</b> to go to the last column.`);
    await this.action('Jump to the first column', this.t!.onCurrentColChanged.pipe(filter(() =>
    this.t!.currentCol.name === columnNames[0])), null, `Press <b>${platformKeyMap['Home'][this.platform]}</b> to return to the first column.`);

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
      this.t!.cell(lastRowIdx - 1, 'RACE').value)), null, `Use <b>${platformKeyMap['Ctrl'][this.platform]}+C</b> to copy the value of the ` +
      `current cell and <b>${platformKeyMap['Ctrl'][this.platform]}+V</b> to paste the copied value into the current cell.`);

    this.title('Selection');

    await this.action(`Select the edited row by clicking on its number (#${lastRowIdx + 1})`,
      this.t!.onSelectionChanged.pipe(filter(() => this.t!.selection.trueCount === 1 && this.t!.selection.get(lastRowIdx))),
      null, 'Grid shares row selection with other viewers, which will be illustrated in the upcoming tutorials. ' +
      'Selected data is highlighted in orange. You can extend selection by holding <b>Shift</b> and dragging the mouse.');

    await this.action('Remove selection by pressing "Esc"', this.t!.onSelectionChanged.pipe(filter(() =>
      !this.t!.selection.anyTrue)), null, 'There are multiple ways to deselect all or certain rows, hitting ' +
      '<b>Esc</b> being the simplest. It removes selection entirely, both from rows and columns.');

    await this.action(`Select all rows with "${platformKeyMap['Ctrl'][this.platform]}+A"`, this.t!.onSelectionChanged.pipe(filter(() =>
      !this.t!.selection.anyFalse)), null, 'There is a key combination to revert this action. You can use ' +
      `<b>${platformKeyMap['Ctrl'][this.platform]}+Shift+A</b> to deselect all rows. It differs from <b>Esc</b>, since it affects only row selection.`);

    const rowCount = this.t!.rowCount;
    const selectedIndexes = [0, 1, 2, 3, 4, rowCount - 5, rowCount - 4, rowCount - 3, rowCount - 2, rowCount - 1];
    await this.action(`Select the first five and the last five table rows (#1 - #5 and #${rowCount - 4} - #${rowCount})`,
      this.t!.onSelectionChanged.pipe(filter(() => this.t!.selection.trueCount === 10 && selectedIndexes.every((idx, i) =>
      this.t!.selection.getSelectedIndexes()[i] === idx))), null, 'Make use of the previously learned shortcuts. ' +
      `<b>Shift+Drag</b> selects rows, <b>${platformKeyMap['Ctrl'][this.platform]}+Shift+Drag</b> deselects rows, <b>${platformKeyMap['Ctrl'][this.platform]}+Click</b> toggles selected state. ` +
      `<b>${platformKeyMap['Ctrl + Shift + Home'][this.platform]}</b> selects rows above the current row, while <b>${platformKeyMap['Ctrl + Shift + End'][this.platform]}</b> selects rows below the ` +
      'current one (inclusively).');

    await this.action('Clear the selection', this.t!.onSelectionChanged.pipe(filter(() => !this.t!.selection.anyTrue)),
      null, `Hit either <b>Esc</b> or <b>${platformKeyMap['Ctrl'][this.platform]}+Shift+A</b>.`);

    grok.shell.windows.showContextPanel = true;
    const noneSeverityRowCount = 3302;
    const severityColumn = this.t!.getCol('SEVERITY');
    await this.action('Find the SEVERITY column and select all rows with the "None" value', this.t!.onSelectionChanged.pipe(
      filter(() => this.t!.selection.trueCount === noneSeverityRowCount)), null, 'First, find any cell with "None" severity ' +
      `and click to make it current. Then press <b>Shift+${platformKeyMap['Enter'][this.platform]}</b>, it will select all rows containing the same value.`);

    await this.action(`Delete these rows (${noneSeverityRowCount})`, this.t!.onRowsRemoved.pipe(filter(() =>
      this.t!.rowCount === (rowCount - noneSeverityRowCount) && !severityColumn.categories.includes('None'))),
      $('div.d4-ribbon-item').has('i.svg-remove-selected-rows')[0], 'In the context panel on the right, you ' +
      'can find the list of actions for selected rows. One of them is "Delete Rows". Pick this option to refine ' +
      'your dataset. You can also use <b>Shift+Delete</b> or the "Remove rows" icon in the top menu.');

    const heightGridCol = grid.col('height')!;
    const weightGridCol = grid.col('weight')!;
    const startedGridCol = grid.col('started')!;
    await this.action('Select HEIGHT and WEIGHT columns', this.t!.onColumnSelectionChanged.pipe(filter(() =>
      heightGridCol.selected && weightGridCol.selected)), null, 'You can do this by holding <b>Shift</b> and clicking ' +
      `the column headers. <b>${platformKeyMap['Ctrl'][this.platform]}</b> will also work in this case. The difference is that holding <b>Shift</b> ` +
      `always adds to selection, whereas <b>${platformKeyMap['Ctrl'][this.platform]}</b> toggles the state (so you can both select and deselect with it).`);

    await this.action('Clear the selection', this.t!.onColumnSelectionChanged.pipe(filter(() => !heightGridCol.selected &&
      !weightGridCol.selected)), null, `Either hit <b>Esc</b> or hold <b>${platformKeyMap['Ctrl'][this.platform]}+Shift</b> while clicking the column headers.`);

    this.title('Sorting, reordering, and resizing');

    await this.action('Sort subjects by AGE in the descending order', grid.onRowsSorted.pipe(filter(() =>
      grid.sortByColumns.length === 1 && grid.sortByColumns[0].name === 'AGE' && grid.sortTypes[0] === false)), null,
      'Double-click the column header once. To change the order to ascending, double-click again. Repeating it one ' +
      'more time will remove the column sorting. You can also sort rows from the column context menu. If a column ' +
      `you want to sort is current, simply press <b>${platformKeyMap['Ctrl'][this.platform]}+Shift+UP</b> in any cell to order its rows.`);

    await this.action('Reset sorting in the grid', interval(1000).pipe(filter(() => grid.sortByColumns.length === 0)),
      null, 'To clear sorting, right-click the column header and select <b>Sort > Reset</b>. The same can be achieved ' +
      'by double-clicking the column header. It switches the state between three modes: descending, ascending, reset.');

    await this.action('Move columns HEIGHT, WEIGHT and STARTED to the beginning of the grid', interval(1000).pipe(
      filter(() => heightGridCol.idx === 1 && weightGridCol.idx === 2 && startedGridCol.idx === 3)), null, 'To reorder ' +
      'columns, either drag each column header to a new position or select the group of columns and move them in one ' +
      'step. When dragging column headers, you can notice two arrows appearing on the grid left and right. Release the ' +
      'mouse over one of them to move the columns to the very beginning of the grid or to the beginning of its visible ' +
      'part. The same is true for the arrows on the right that help you move columns to the grid end.<br>Make sure ' +
      'that the new order is HEIGHT, WEIGHT, STARTED, followed by all other columns.');

    const initialRowHeight = grid.props.rowHeight;
    await this.action('Increase row height', interval(1000).pipe(filter(() => grid.props.rowHeight > initialRowHeight)),
      null, 'Dragging a row header border will resize the row height (move up to decrease it, or down to increase it).');

    const initialColWidth = startedGridCol.width;
    await this.action('Extend the STARTED column width', interval(1000).pipe(filter(() => startedGridCol.width >
      initialColWidth)), null, 'To resize a column width, drag the column header. This also works for multiple selected ' +
      'columns. You can also adjust column sizing for the whole grid: right-click on any cell in the grid and select ' +
      '<b>Column Sizing</b> (options are Minimal, Compact, Optimal, or Maximal).');

    this.title('Formatting');

    const heightColumn = this.t!.getCol('height');
    const weightColumn = this.t!.getCol('weight');
    const startedColumn = this.t!.getCol('started');
    await this.action('Change number formatting of HEIGHT to scientific notation', this.t!.onMetadataChanged.pipe(
      filter((data) => (data.args.source as unknown as DG.Column).name === heightColumn.name &&
      (data.args.key as unknown as string) === 'format' && (data.args.change as unknown as string) === 'set' &&
      (data.args.value as unknown as string) === 'scientific')), null, 'One way to set a new numeric format is to ' +
      'right-click a column header and select <b>Format > scientific</b>. For convenience, there are examples for ' +
      'each format option in the menu. Hover the mouse over the list of options to see the format name shown in a tooltip.');

    await this.action('Set number formatting of WEIGHT to "max two digits after comma"', this.t!.onMetadataChanged.pipe(
      filter((data) => (data.args.source as unknown as DG.Column).name === weightColumn.name && (data.args.key as unknown as
      string) === 'format' && (data.args.change as unknown as string) === 'set' && (data.args.value as unknown as string)
      === 'max two digits after comma')));

    await this.action('Set date formatting of STARTED to "dd.MM.yyyy"', this.t!.onMetadataChanged.pipe(filter((data) =>
      (data.args.source as unknown as DG.Column).name === startedColumn.name && (data.args.key as unknown as string) ===
      'format' && (data.args.change as unknown as string) === 'set' && (data.args.value as unknown as string) ===
      'dd.MM.yyyy')), null, 'Learn more about ' + ui.link('formatting of numbers and dates',
      'https://datagrok.ai/help/discover/tags#format').outerHTML);

    this.title('Color coding');

    await this.action('Set linear color coding for HEIGHT', this.t!.onMetadataChanged.pipe(filter((data) =>
      (data.args.source as unknown as DG.Column).name === heightColumn.name && (data.args.key as unknown as string) ===
      '.color-coding-type' && (data.args.change as unknown as string) === 'set' && (data.args.value as unknown as string) ===
      'Linear')), null, 'Right-click the column header and enable <b>Color Coding > Linear</b>. This action will paint the ' +
      'column cells in colors ranging from blue to red (by default). This way you can quickly tell which values are closer ' +
      'to the column minimum and which are closer to the maximum. There are other color coding types (categorical and ' +
      'conditional) as well as various color coding related options. Learn more about this topic in ' +
      ui.link('the wiki article', 'https://datagrok.ai/help/visualize/viewers/grid#color-coding').outerHTML);

    this.title('Summary columns');

    await this.action('Add a summary column with inline bar chart', interval(1000).pipe(filter(() => grid.col('barchart') !==
      null)), null, 'Summary columns is a way to visualize multiple numerical values across the row. To add one, right-click ' +
      'on any cell in the grid and select <b>Add > Summary Columns > Bar Chart</b>.');

    await this.action('Click the "barchart" column header', interval(1000).pipe(filter(() => grok.shell.o instanceof
      DG.GridColumn && grok.shell.o.name === 'barchart')), null, 'Summary columns contain an inline viewer that visualizes by ' +
      'default up to 3 columns. When you click the column header, you can see which columns are used for a summary in the ' +
      'context panel on the right.');
  }
}
