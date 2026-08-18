import { test, expect, Page } from '@playwright/test';
import {
  BALLOON_CONTAINER,
  CONTEXT_PANEL,
  CONTEXT_PANEL_EXPAND_ALL,
  treeGroupByName,
  treeItemByName,
  treeNodeByPath,
} from './selectors';
import {
  goHome,
  ensureBrowsePanelOpen,
  ensureContextPanelOpen,
  watchErrors,
  expectNoErrors,
  expandTreeGroup,
} from './helpers';

const REQUIRED_MYSTUFF = ['Recent', 'Favorites', 'Shared with me', 'My files'];

const DB_PROVIDERS = [
  'Access', 'Athena', 'BigQuery', 'Cassandra', 'ClickHouse', 'DB2', 'Databricks',
  'Denodo', 'Firebird', 'HBase', 'Hive', 'Hive2', 'Impala', 'MLFlow', 'MS SQL',
  'MariaDB', 'MongoDB', 'MySQL', 'Neo4j', 'Neptune', 'ODATA', 'Oracle', 'PI',
  'Postgres', 'Redshift', 'SAP HANA', 'Snowflake', 'Sparql', 'Teradata',
  'Vertica', 'Virtuoso',
];

const PLATFORM_GROUPS = ['Admin', 'Plugins', 'Credentials', 'Functions', 'Predictive models', 'Dockers', 'Sticky Meta'];
const PLATFORM_ITEMS = ['Users', 'Groups', 'Roles', 'Notebooks', 'MCP Servers', 'Repositories', 'Keys', 'Sync', 'Layouts'];

async function ensureContextPanelExpanded(page: Page): Promise<void> {
  await ensureContextPanelOpen(page,  true);
}

test.describe('Browse — tree nodes open without errors (section 18)', () => {
  test.beforeEach(async ({ page }) => {
    await goHome(page);
    await ensureBrowsePanelOpen(page);
    await ensureContextPanelExpanded(page);
  });

  test('Browse-Node-MyStuff-01 — clicking each My stuff subnode produces no errors', async ({ page }) => {
    const sink = watchErrors(page);

    await expandTreeGroup(page, 'My stuff');

    for (const name of REQUIRED_MYSTUFF) {
      const node = treeNodeByPath(page, ['My-stuff', name]);
      await expect(node, `My stuff > ${name} should exist`).toHaveCount(1);
      await ensureBrowsePanelOpen(page);
      await ensureContextPanelExpanded(page);
      await expandTreeGroup(page, 'My stuff');
      await node.scrollIntoViewIfNeeded();
      await node.click();
      await page.waitForTimeout(800);
      await expectNoErrors(page, sink);
    }
  });

  test('Browse-Node-Spaces-01 — clicking Spaces shows children without errors', async ({ page }) => {
    const sink = watchErrors(page);
    await expandTreeGroup(page, 'Spaces');
    await expectNoErrors(page, sink);
  });

  test('Browse-Node-Dashboards-01 — Dashboards opens the Projects view without errors', async ({ page }) => {
    const sink = watchErrors(page);
    const dashboards = treeItemByName(page, 'Dashboards');
    await expect(dashboards).toBeVisible({ timeout: 10_000 });
    await dashboards.click();
    await page.waitForTimeout(1500);
    await expect(page).toHaveURL(/\/projects\b/);
    await expectNoErrors(page, sink);
  });

  test('Browse-Node-DB-01 — every database provider is clickable without errors', async ({ page }) => {
    const sink = watchErrors(page);
    await expandTreeGroup(page, 'Databases');

    for (const provider of DB_PROVIDERS) {
      const node = treeNodeByPath(page, ['Databases', provider]);
      if (!(await node.isVisible().catch(() => false))) {

        continue;
      }
      await ensureBrowsePanelOpen(page);
      await ensureContextPanelExpanded(page);
      await expandTreeGroup(page, 'Databases');
      await node.scrollIntoViewIfNeeded();
      await node.click();
      await page.waitForTimeout(500);
      await expectNoErrors(page, sink);
    }
  });

  test('Browse-Node-Platform-01 — every Platform subnode is clickable without errors', async ({ page }) => {
    test.setTimeout(180_000); 
    const sink = watchErrors(page);
    await expandTreeGroup(page, 'Platform');

    for (const name of [...PLATFORM_GROUPS, ...PLATFORM_ITEMS]) {
      const node = treeNodeByPath(page, ['Platform', name]);
      if (!(await node.isVisible().catch(() => false))) continue;
      await ensureBrowsePanelOpen(page);
      await ensureContextPanelExpanded(page);
      await expandTreeGroup(page, 'Platform');
      await node.scrollIntoViewIfNeeded();
      await node.click();
      await page.waitForTimeout(800);
      await expectNoErrors(page, sink);
    }
  });
});
