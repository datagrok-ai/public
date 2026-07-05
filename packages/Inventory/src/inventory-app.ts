import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {inventoryDb, ItemsRow} from './generated/db';

const reasons = ['received', 'shipped', 'adjustment', 'damaged', 'returned'];
const systemColumns = ['id', 'version', 'created_on', 'updated_on', 'author_id'];

/** Adjusts the stock of [itemId] by [delta] using optimistic concurrency: reads the
 * fresh row, then atomically updates `quantity` (guarded by its `version`) and logs
 * the `stock_movements` row in one transaction — either both land or neither does.
 * Retries with a fresh read on a version conflict (a concurrent adjustment bumped
 * the version between the read and the write). */
export async function adjustStock(itemId: string, delta: number, reason: string,
  maxRetries: number = 5): Promise<ItemsRow> {
  for (let attempt = 0; ; attempt++) {
    const item = await inventoryDb.items.get(itemId);
    if (item == null)
      throw new Error('Item not found or not visible');
    const quantity = (item.quantity ?? 0) + delta;
    try {
      const [updated] = await grok.dapi.domains.transaction('inventory', [
        {op: 'update', table: 'items', id: itemId, values: {quantity: quantity}, expectedVersion: item.version},
        {op: 'insert', table: 'stock_movements',
          values: {item_id: itemId, delta: delta, reason: reason, moved_on: new Date().toISOString()}},
      ]);
      return {...item, quantity: quantity, version: updated.version};
    } catch (e: any) {
      if (attempt >= maxRetries || !`${e?.message ?? e}`.includes('Version conflict'))
        throw e;
    }
  }
}

/** Minimal inventory app over the typed clients that `grok api` generates from
 * databases/inventory/schema.json (src/generated/db.ts): items grid via queryDf
 * (d42 column tags drive the property panel and renderers), CSV/Parquet import
 * with batch upsert by SKU, optimistic stock adjustments, and a
 * stock-on-hand-by-location aggregation view. */
export class InventoryApp {
  view: DG.ViewBase = DG.View.create();
  itemsDf: DG.DataFrame | null = null;
  gridHost = ui.div([]);
  aggHost = ui.div([]);
  private refreshSeq = 0;

  static async run(): Promise<DG.ViewBase> {
    const app = new InventoryApp();
    await app.init();
    return app.view;
  }

  async init(): Promise<void> {
    this.view.name = 'Inventory';
    this.view.setRibbonPanels([[
      ui.bigButton('Import...', () => this.importDialog(), 'Import a CSV or Parquet file; rows are upserted by SKU'),
      ui.button('Adjust stock...', () => this.adjustDialog()),
      ui.iconFA('sync', () => this.refresh(), 'Refresh'),
    ]]);

    const tabs = ui.tabControl();
    tabs.addPane('Items', () => this.gridHost);
    tabs.addPane('Stock by location', () => {
      this.refreshAggregation();
      return this.aggHost;
    });
    tabs.onTabChanged.subscribe(() => {
      if (tabs.currentPane?.name === 'Stock by location')
        this.refreshAggregation();
    });
    tabs.root.style.width = '100%';
    tabs.root.style.height = '100%';
    this.view.root.appendChild(tabs.root);
    await this.refresh();
  }

  async refresh(): Promise<void> {
    const seq = ++this.refreshSeq;
    let df: DG.DataFrame;
    try {
      df = await inventoryDb.items.queryDf({sort: 'sku'});
    } catch (e: any) {
      grok.shell.error(e.message ?? `${e}`);
      return;
    }
    if (seq !== this.refreshSeq)
      return;
    this.itemsDf = df;
    ui.empty(this.gridHost);
    if (df.rowCount === 0) {
      this.gridHost.appendChild(ui.divText(
        'No items yet — import a CSV or Parquet file to get started.'));
      return;
    }
    df.name = 'Inventory items';
    const grid = DG.Viewer.grid(df);
    for (const c of systemColumns)
      if (grid.col(c) != null)
        grid.col(c)!.visible = false;
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.gridHost.appendChild(grid.root);
  }

  async refreshAggregation(): Promise<void> {
    ui.empty(this.aggHost);
    this.aggHost.appendChild(ui.loader());
    let rows: any[];
    try {
      rows = await inventoryDb.items.aggregate({
        groupBy: ['location'],
        measures: [
          {fn: 'sum', column: 'quantity', as: 'stock_on_hand'},
          {fn: 'count', as: 'items'},
        ],
        sort: 'location',
      });
    } catch (e: any) {
      ui.empty(this.aggHost);
      this.aggHost.appendChild(ui.divText(e.message ?? `${e}`));
      return;
    }
    ui.empty(this.aggHost);
    if (rows.length === 0) {
      this.aggHost.appendChild(ui.divText('No items to aggregate.'));
      return;
    }
    const df = DG.DataFrame.fromObjects(rows)!;
    df.name = 'Stock on hand by location';
    const grid = DG.Viewer.grid(df);
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    this.aggHost.appendChild(grid.root);
  }

  importDialog(): void {
    DG.Utils.openFile({
      accept: '.csv,.parquet',
      open: (file) => this.confirmImport(file),
    });
  }

  confirmImport(file: File): void {
    const partial = ui.input.bool('Apply good rows', {value: false});
    partial.setTooltip('Apply valid rows and report the bad ones per row instead of aborting the whole file');
    ui.dialog('Import items')
      .add(ui.divV([
        ui.divText(`${file.name} → inventory.items, upsert by SKU`),
        ui.inputs([partial]),
      ]))
      .onOK(() => this.runImport(file, partial.value === true))
      .show();
  }

  async runImport(file: File, partial: boolean): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create(`Importing ${file.name}...`);
    try {
      const options = {mode: 'upsert' as const, allOrNothing: !partial};
      const report = file.name.toLowerCase().endsWith('.parquet')
        ? await inventoryDb.items.batch(new Uint8Array(await file.arrayBuffer()),
          {...options, format: 'parquet'})
        : await inventoryDb.items.batch(await file.text(), options);
      const summary = `${report.inserted} inserted, ${report.updated} updated, ${report.skipped} skipped`;
      if (report.error != null || report.errorCount > 0)
        grok.shell.warning(`Import ${report.error != null ? 'aborted' : 'finished with errors'}: ` +
          `${report.errorCount} bad rows; ${summary}`);
      else
        grok.shell.info(`Import complete: ${summary}`);
      await this.refresh();
    } catch (e: any) {
      grok.shell.error(e.message ?? `${e}`);
    } finally {
      pi.close();
    }
  }

  adjustDialog(): void {
    const item = this.currentItem();
    if (item == null) {
      grok.shell.warning('Select an item in the grid first');
      return;
    }
    const delta = ui.input.int('Delta', {value: 0});
    delta.setTooltip('Signed quantity change; the result may not go below zero');
    const reason = ui.input.choice('Reason', {items: reasons, value: 'adjustment'});
    ui.dialog(`Adjust stock: ${item.sku}`)
      .add(ui.inputs([delta, reason]))
      .onOK(async () => {
        if (!delta.value) {
          grok.shell.warning('Delta must be a non-zero number');
          return;
        }
        try {
          const updated = await adjustStock(item.id, delta.value, reason.value ?? 'adjustment');
          grok.shell.info(`${item.sku}: ${updated.quantity} on hand`);
          await this.refresh();
        } catch (e: any) {
          grok.shell.error(e.message ?? `${e}`);
        }
      })
      .show();
  }

  currentItem(): {id: string, sku: string} | null {
    const df = this.itemsDf;
    if (df == null || df.rowCount === 0)
      return null;
    const idx = df.currentRowIdx >= 0 ? df.currentRowIdx : 0;
    return {id: df.col('id')!.get(idx), sku: df.col('sku')!.get(idx)};
  }
}
