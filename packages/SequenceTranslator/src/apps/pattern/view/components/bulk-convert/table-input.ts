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

  private getTableFromEventData(eventData: any): DG.DataFrame {
    if (!eventData && eventData.args && eventData.args.dataFrame instanceof DG.DataFrame)
      throw new Error(`EventData does not contain a dataframe`, eventData);

    return eventData.args.dataFrame as DG.DataFrame;
  }

  private handleTableAdded(eventData: any): void {
    const table = this.getTableFromEventData(eventData);

    if (this.availableTables.some((t: DG.DataFrame) => t.name === table.name))
      return;

    this.availableTables.push(table);
    this.eventBus.selectTable(table);
    this.refreshTableInput();
  }

  private handleTableRemoved(eventData: any): void {
    const removedTable = this.getTableFromEventData(eventData);
    this.availableTables = this.availableTables.filter((table: DG.DataFrame) => table.name !== removedTable.name);

    const table = this.availableTables[0];
    this.eventBus.selectTable(table ? table : null);
    this.refreshTableInput();
  }

  private refreshTableInput(): void {
    const tableInput = this.createTableInput();
    $(this.tableInputContainer).empty();
    this.tableInputContainer.append(tableInput.root);
  }

  private createTableInput(): DG.InputBase<DG.DataFrame | null> {
    const currentlySelectedTable = this.eventBus.getTableSelection();

    const tableInput = ui.input.table('Tables', {
      value: currentlySelectedTable!, items: this.availableTables,
      onValueChanged: (value) => {
        // WARNING: non-null check necessary to prevent resetting columns to
        // null upon handling onTableAdded
        if (value !== null)
          this.eventBus.selectTable(value);
      }
    });
    return tableInput;
  }

  private handleTableChoice(): void {
    const selectedTable = this.eventBus.getTableSelection();
    if (!selectedTable) return;
    if (!this.isTableDisplayed(selectedTable))
      this.displayTable(selectedTable);
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
