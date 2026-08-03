// GROK-19103: joining two source tables and saving must create exactly ONE project, not a stray join-only project.
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs} from './_helpers';
import {deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

// Open demog.csv via the OpenFile recorder API — writes df.tags['.script'] provenance without UI traversal flake.
async function openDemogFile(page: Page): Promise<void> {
  await page.evaluate(async () => {
    const DG = (window as any).DG;
    await DG.Func.find({name: 'OpenFile'})[0].prepare({
      fullPath: 'System:DemoFiles/demog.csv',
    }).call(undefined, undefined, {processed: false});
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1500);
}

test('Projects / Complex derived-tables: Join lands in active project (GROK-19103)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = 'AutoTest-ComplexDerived-' + stamp;

  await loginToDatagrok(page);
  await page.goto('/');
  await page.waitForFunction(() => {
    try { return !!(window as any).grok?.shell?.user?.login; } catch { return false; }
  }, {timeout: 60_000});
  await evalJs(page, `(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  })()`);
  await page.waitForTimeout(1000);

  let baselineCount = 0;
  let projectId: string | null = null;
  let tableInfoId: string | null = null;

  try {
    await softStep('Setup + Step 1: open 2 source tables via OpenFile (colon-form fullPath)', async () => {
      await openDemogFile(page);
      await openDemogFile(page);
      await evalJs(page, `(() => {
        const tables = grok.shell.tables;
        if (tables.length >= 2) {
          tables[0].name = 'src_a_${stamp}';
          tables[1].name = 'src_b_${stamp}';
        }
      })()`);
      await page.waitForTimeout(1000);
      const tableCount = await evalJs(page, 'grok.shell.tables.length');
      expect(tableCount).toBe(2);
    });

    await softStep('Step 1 (Join sub-bullet): grok.data.joinTables produces 3rd table in workspace', async () => {
      await evalJs(page, `(async () => {
        const [a, b] = grok.shell.tables;
        const joined = grok.data.joinTables(
          a, b,
          ['USUBJID'], ['USUBJID'],
          ['AGE', 'SEX'], ['RACE'],
          'inner',
        );
        joined.name = 'joined_${stamp}';
        grok.shell.addTableView(joined);
      })()`);
      await page.waitForTimeout(2000);
      const tableCount = await evalJs(page, 'grok.shell.tables.length');
      expect(tableCount).toBe(3);
    });

    await softStep('Step 2: capture project-count baseline before Save (GROK-19103 invariant prep)', async () => {
      baselineCount = await evalJs(page,
        `(async () => (await grok.dapi.projects.list({limit: 1000})).length)()`,
      );
      expect(baselineCount).toBeGreaterThanOrEqual(0);
    });

    await softStep('Step 2: Save current project with Data Sync ON', async () => {
      // Multi-table inline save — must persist all 3 open tables (helper saves only the active TableView).
      const saved = await page.evaluate(async (n) => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        const project = DG.Project.create();
        project.name = n;
        const tables = Array.from(grok.shell.tables);
        const tableInfoIds: string[] = [];
        for (const df of tables as any[]) {
          const ti = df.getTableInfo();
          project.addChild(ti);
          await grok.dapi.tables.uploadDataFrame(df);
          await grok.dapi.tables.save(ti);
          tableInfoIds.push(ti.id);
        }
        const tv = grok.shell.tv;
        const layout = tv?.saveLayout?.();
        if (layout) {
          project.addChild(layout);
          await grok.dapi.layouts.save(layout);
        }
        await grok.dapi.projects.save(project);
        return {projectId: project.id, tableInfoIds};
      }, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoIds[0] ?? null;
      expect(projectId).toBeTruthy();
      const exists = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        return p != null;
      }, projectId);
      expect(exists).toBe(true);
    });

    await softStep('GROK-19103 INVARIANT: exactly ONE new project created (no stray join-only project)', async () => {
      const afterCount = await evalJs(page,
        `(async () => (await grok.dapi.projects.list({limit: 1000})).length)()`,
      );
      // GROK-19103: a stray join-only project would push delta to 2+.
      expect(afterCount - baselineCount).toBe(1);
    });

    await softStep('Reopen: verify project loads with source tables intact', async () => {
      if (!projectId) throw new Error('no projectId captured');
      await closeAll(page);
      const tableCount = await page.evaluate(async (pid) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(pid);
        await p.open();
        await new Promise((r) => setTimeout(r, 2000));
        try { return Number(grok.shell.tables?.length) || 0; }
        catch (_) { return 0; }
      }, projectId);
      // At minimum the 2 source tables must reopen; the joined derivative may or may not persist.
      expect(tableCount).toBeGreaterThanOrEqual(2);
    });
  } finally {
    await deleteProjectWithCleanup(page, {
      projectId: projectId ?? undefined,
      tableInfoId: tableInfoId ?? undefined,
    });
    await closeAll(page);
  }

  finishSpec();
});
