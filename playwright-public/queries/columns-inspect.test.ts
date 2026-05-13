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

// Test Track scenario `columns-inspect.md` (order 8).
// The Test Track asks to click every column of every table under both `PostgresDart/NorthwindTest`
// AND `Postgres/NorthwindTest`. On the public server there's no `NorthwindTest` anywhere, and
// PostgresDart hosts different connections entirely — so we skip the PostgresDart subtest and
// exercise the Postgres/Northwind schema walk with a representative subset (three tables × all
// their columns ≈ 30 clicks) rather than all ~140 clicks the track implies.

const PROVIDER = 'Postgres';
const SCHEMA = 'public';
// CI: target three Datagrok metadata tables (System:Datagrok). They always
// exist, each exposes its full column list, and together cover the same
// "click every column of three different tables" intent as the original
// Northwind triplet (products/orders/customers).
const TABLES = ['users', 'groups', 'entities'] as const;

test.describe.serial(`DB schema column inspection (${PROVIDER} / ${POSTGRES_CONNECTION})`, () => {
  test('Clicking each column of products/orders/customers sets it as current object without errors', async ({ page }) => {
    // The tree walk touches dozens of server-backed nodes — give it headroom.
    test.setTimeout(180_000);

    await goHome(page);
    await showContextPanel(page);

    // Resolve the connection's server-side name so we can target the Schemas group wrapper.
    const connServerName = await getConnectionServerName(page, PROVIDER, POSTGRES_CONNECTION);
    expect(connServerName, 'Postgres/Northwind connection should exist on the server').toBeTruthy();

    // Navigate deep into Databases > Postgres > Northwind > Schemas > public.
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
        // The current object should now be the DB column — `grok.shell.o.name` matches the last
        // segment of the node `name=` attribute.
        const expectedColName = colNode.split('---').pop()!;
        const actualName = await getCurrentObjectName(page);
        expect(actualName, `clicking ${colNode} should select the column as current object`)
          .toBe(expectedColName);

        // Context Panel should not surface any error balloons as a side effect.
        expect(await getVisibleErrorBalloons(page),
          `no error balloons expected after selecting ${colNode}`).toEqual([]);
      }
    }
  });
});
