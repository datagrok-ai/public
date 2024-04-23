/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import '../../style.css';

import {EventBus} from '../../../model/event-bus';

export class TableInputManager {
  private availableTables: DG.DataFrame[] = [];
  private tableInputContainer: HTMLDivElement = ui.div([]);

  constructor(private eventBus: EventBus) {
    this.subscribeToTableEvents();
    this.refreshTableInput();
  }

  getTableInputContainer(): HTMLDivElement {
    return this.tableInputContainer;
  }

  private subscribeToTableEvents(): void {
    grok.events.onTableAdded.subscribe((eventData) => this.handleTableAdded(eventData));
    grok.events.onTableRemoved.subscribe((eventData) => this.handleTableRemoved(eventData));
    this.eventBus.tableSelectionChanged$.subscribe(() => this.handleTableChoice());
  }

  private handleTableAdded(eventData: DG.EventData<DG.DataFrameArgs>): void {
    const table = eventData.args.dataFrame;
    if (this.availableTables.some((availableTable: DG.DataFrame) => availableTable.name === table.name))
      return;

    this.availableTables.push(table);

    this.refreshTableInput();
  }

  private handleTableRemoved(eventData: DG.EventData<DG.DataFrameArgs>): void {
    const removedTable = eventData.args.dataFrame;
    this.availableTables = this.availableTables.filter((table: DG.DataFrame) => table.name !== removedTable.name);

    this.refreshTableInput();
  }

  private refreshTableInput(): void {
    const tableInput = this.createTableInput();
    $(this.tableInputContainer).empty();
    $(this.tableInputContainer).append(tableInput.root);
  }

  private createTableInput(): DG.InputBase<DG.DataFrame | null> {
    const currentSelection = this.eventBus.getTableSelection() ?
      this.eventBus.getTableSelection() : this.availableTables[0];

    const tableInput = ui.tableInput(
      'Tables',
      currentSelection,
      this.availableTables
    );
    tableInput.onChanged(
      () => this.eventBus.selectTable(tableInput.value)
    );

    this.eventBus.selectTable(currentSelection);

    return tableInput;
  }

  private handleTableChoice(): void {
    const table = this.eventBus.getTableSelection();
    if (!table) return;
    if (!this.isTableDisplayed(table))
      this.displayTable(table);
  }

  private isTableDisplayed(table: DG.DataFrame): boolean {
    return grok.shell.tableNames.includes(table.name);
  }

  private displayTable(table: DG.DataFrame): void {
    const previousView = grok.shell.v;
    grok.shell.addTableView(table);
    grok.shell.v = previousView;
  }
}
