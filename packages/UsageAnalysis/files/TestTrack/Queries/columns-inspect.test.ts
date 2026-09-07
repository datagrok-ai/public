import { test, expect } from '@playwright/test';
import {
  POSTGRES_CONNECTION,
  expandDbConnection,
  expandDbProvider,
  expandDbSchemas,
  expandTreeNode,
  getConnectionServerName,
  getCurrentObjectName,
  getVisibleErrorBalloons,
  goHome,
  listDbTableColumnNodeNames,
  selectTreeNodeAsCurrentObject,
  showContextPanel,
} from './helpers';

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
const TABLES = ['products', 'orders', 'customers'] as const;

test.describe.serial(`DB schema column inspection (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test('Clicking each column of products/orders/customers sets it as current object without errors', async ({ page }) => {

    test.setTimeout(180_000);

    await goHome(page);
    await showContextPanel(page);

    const connServerName = await getConnectionServerName(page, PROVIDER, POSTGRES_CONNECTION);
    expect(connServerName, 'Postgres/Northwind connection should exist on the server').toBeTruthy();

    await expandDbProvider(page, PROVIDER);
    await expandDbConnection(page, PROVIDER, POSTGRES_CONNECTION);
    await expandDbSchemas(page, PROVIDER, connServerName);
    const schemaNode =
      `tree-Databases---${PROVIDER}---${POSTGRES_CONNECTION}---Schemas---${SCHEMA}`;
    await expandTreeNode(page, schemaNode);

    for (const table of TABLES) {
      const tableNode = `${schemaNode}---${table}`;
      await expandTreeNode(page, tableNode);

      const columnNodes = await listDbTableColumnNodeNames(page, tableNode);
      expect(columnNodes.length, `table ${table} should expose columns`).toBeGreaterThan(0);

      for (const colNode of columnNodes) {
        await selectTreeNodeAsCurrentObject(page, colNode);

        const expectedColName = colNode.split('---').pop()!;
        const actualName = await getCurrentObjectName(page);
        expect(actualName, `clicking ${colNode} should select the column as current object`)
          .toBe(expectedColName);

        expect(await getVisibleErrorBalloons(page),
          `no error balloons expected after selecting ${colNode}`).toEqual([]);
      }
    }
  });
});
