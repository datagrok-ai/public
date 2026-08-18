import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Scripts Layout — test_Layout: Layout pane viewers + project round-trip", async ({
  page,
}) => {
  test.setTimeout(360_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add("selenium");
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.windows.showToolbox = true;
    grok.shell.closeAll();
    const prior = await grok.dapi.scripts.filter('name = "test_Layout"').list();
    for (const s of prior)
      try {
        await grok.dapi.scripts.delete(s);
      } catch (_) {}
    const priorProj = await grok.dapi.projects
      .filter('name = "test_Layout_proj_v3"')
      .list();
    for (const p of priorProj)
      try {
        await grok.dapi.projects.delete(p);
      } catch (_) {}
  });

  let scriptId: string | null = null;
  let projectId: string | null = null;
  let layoutId: string | null = null;

  await softStep("1. Create test_Layout JavaScript script", async () => {
    scriptId = await page.evaluate(async () => {
      const body = `//name: test_Layout
//language: javascript
//input: int idx=1
//output: dataframe df

const fileList = await grok.dapi.files.list('System:DemoFiles/chem', true, '');
const csvFiles = fileList.filter((fi) => fi.fileName.endsWith('.csv'));
if (csvFiles.length === 0)
  throw new Error('No CSV files found in System:DemoFiles/chem');

csvFiles.sort((a, b) => b.createdOn - a.createdOn);
const lastModifiedFile = csvFiles[idx];
const csv = await grok.dapi.files.readAsText(lastModifiedFile.fullPath);
df = DG.DataFrame.fromCsv(csv);
`;
      const s = DG.Script.create(body);
      const saved = await grok.dapi.scripts.save(s);
      return saved.id;
    });
    expect(scriptId).toBeTruthy();
  });

  await softStep("2. Open script editor, switch to Layout tab", async () => {
    await page.evaluate((id) => {
      grok.shell.route(`/script/${id}`);
    }, scriptId);
    await page.waitForFunction(
      () => (window as any).grok?.shell?.v?.name === "test_Layout",
      null,
      { timeout: 30000 },
    );
    await page.evaluate(async () => {
      const tab = Array.from(document.querySelectorAll(".d4-tab-header")).find(
        (e) =>
          e.textContent?.trim() === "Layout" &&
          (e as HTMLElement).offsetParent !== null,
      ) as HTMLElement;
      tab?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    const hasRunLink = await page.evaluate(() =>
      Array.from(document.querySelectorAll("a, span, .d4-link-label, button"))
        .filter((e) => (e as HTMLElement).offsetParent !== null)
        .some((e) => e.textContent?.trim() === "Run script"),
    );
    expect(hasRunLink).toBe(true);
  });

  await softStep(
    "3. Click Run script; OK in dialog; grid appears",
    async () => {
      await page.evaluate(async () => {
        const runLink = Array.from(
          document.querySelectorAll("a, span, .d4-link-label, button"),
        )
          .filter((e) => (e as HTMLElement).offsetParent !== null)
          .find((e) => e.textContent?.trim() === "Run script") as HTMLElement;
        runLink?.click();
        for (let i = 0; i < 60; i++) {
          const ok = Array.from(
            document.querySelectorAll('[name="button-OK"]'),
          ).find(
            (b) => (b as HTMLElement).offsetParent !== null,
          ) as HTMLElement;
          if (ok) {
            ok.click();
            break;
          }
          await new Promise((r) => setTimeout(r, 100));
        }
        for (let i = 0; i < 100; i++) {
          if (document.querySelector('[name="viewer-Grid"]')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 1500));
      });
      const hasGrid = await page.evaluate(
        () => !!document.querySelector('[name="viewer-Grid"]'),
      );
      expect(hasGrid).toBe(true);
    },
  );

  await softStep(
    "4. Add Scatter Plot + Histogram via Toolbox viewer icons (Layout pane)",
    async () => {

      const viewerNames = await page.evaluate(async () => {
        grok.shell.windows.simpleMode = false;
        grok.shell.windows.showToolbox = true;
        await new Promise((r) => setTimeout(r, 800));

        const toolboxSpan = Array.from(
          document.querySelectorAll(".tab-handle-text"),
        ).find((e) => e.textContent?.trim() === "Toolbox") as
          | HTMLElement
          | undefined;
        (toolboxSpan?.parentElement as HTMLElement)?.click();
        await new Promise((r) => setTimeout(r, 800));
        (
          document.querySelector('[name="icon-scatter-plot"]') as HTMLElement
        )?.click();
        await new Promise((r) => setTimeout(r, 1500));
        (
          document.querySelector('[name="icon-histogram"]') as HTMLElement
        )?.click();
        await new Promise((r) => setTimeout(r, 1500));
        return Array.from(document.querySelectorAll('[name^="viewer-"]')).map(
          (e) => e.getAttribute("name"),
        );
      });
      expect(viewerNames).toContain("viewer-Grid");
      expect(viewerNames).toContain("viewer-Scatter-plot");
      expect(viewerNames).toContain("viewer-Histogram");
    },
  );

  await softStep(
    "6. Save script — saveScript walks layoutViews and persists the layout",
    async () => {

      const layoutOnParam = await page.evaluate(async () => {
        const save = Array.from(
          document.querySelectorAll('[name="button-Save"]'),
        ).find((b) => (b as HTMLElement).offsetParent !== null) as HTMLElement;
        save?.click();
        for (let i = 0; i < 30; i++) {
          await new Promise((r) => setTimeout(r, 500));
          const s = await grok.dapi.scripts
            .filter('name = "test_Layout"')
            .first();
          const id = s?.outputs?.[0]?.options?.["layout"];
          if (id) return id;
        }
        return null;
      });
      expect(layoutOnParam).toBeTruthy();
    },
  );

  await softStep("7. Close All", async () => {
    await page.evaluate(async () => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1000));
    });
    const view = await page.evaluate(() => grok.shell.v?.name);
    expect(view).toBe("Home");
  });

  await softStep(
    "8. Re-run via editor Layout tab — saved layout reapplies (Grid+SP+Hist)",
    async () => {

      const viewerNames = await page.evaluate(async (id) => {
        grok.shell.route(`/script/${id}`);
        for (let i = 0; i < 60; i++) {
          if (grok.shell.v?.name === "test_Layout") break;
          await new Promise((r) => setTimeout(r, 500));
        }
        await new Promise((r) => setTimeout(r, 1000));
        const tab = Array.from(
          document.querySelectorAll(".d4-tab-header"),
        ).find(
          (e) =>
            e.textContent?.trim() === "Layout" &&
            (e as HTMLElement).offsetParent !== null,
        ) as HTMLElement;
        tab?.click();
        await new Promise((r) => setTimeout(r, 1500));
        const runLink = Array.from(
          document.querySelectorAll("a, span, .d4-link-label, button"),
        )
          .filter((e) => (e as HTMLElement).offsetParent !== null)
          .find((e) => e.textContent?.trim() === "Run script") as HTMLElement;
        runLink?.click();
        for (let i = 0; i < 60; i++) {
          const ok = Array.from(
            document.querySelectorAll('[name="button-OK"]'),
          ).find(
            (b) => (b as HTMLElement).offsetParent !== null,
          ) as HTMLElement;
          if (ok) {
            ok.click();
            break;
          }
          await new Promise((r) => setTimeout(r, 100));
        }
        for (let i = 0; i < 100; i++) {
          if (document.querySelector('[name="viewer-Grid"]')) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 4000));
        return Array.from(document.querySelectorAll('[name^="viewer-"]')).map(
          (e) => e.getAttribute("name"),
        );
      }, scriptId);
      expect(viewerNames).toContain("viewer-Grid");
      expect(viewerNames).toContain("viewer-Scatter-plot");
      expect(viewerNames).toContain("viewer-Histogram");
    },
  );

  await softStep(
    "9. Add new viewers via JS API on script-output TableView",
    async () => {

      const types = await page.evaluate(async () => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1000));
        const s = await grok.dapi.scripts
          .filter('name = "test_Layout"')
          .first();
        const res = await s.apply({ idx: 1 });
        res.name = "df_layout";
        grok.shell.addTableView(res);
        await new Promise((r) => setTimeout(r, 2000));
        const tv = grok.shell.tv;
        tv.addViewer(DG.VIEWER.SCATTER_PLOT);
        tv.addViewer(DG.VIEWER.HISTOGRAM);
        await new Promise((r) => setTimeout(r, 1500));
        return Array.from(tv.viewers as any).map((v: any) => v?.type);
      });
      expect(types).toContain("Grid");
      expect(types).toContain("Scatter plot");
      expect(types).toContain("Histogram");
    },
  );

  await softStep(
    "10. Save the project (dataframe + linked layout)",
    async () => {
      const ids = await page.evaluate(async () => {
        const tv = grok.shell.tv;
        const project = DG.Project.create();
        project.name = "test_Layout_proj_v3";
        const tableInfo = tv.dataFrame.getTableInfo();
        project.addChild(tableInfo);
        const layout = tv.saveLayout();
        await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
        await grok.dapi.tables.save(tableInfo);
        await grok.dapi.layouts.save(layout);
        project.addChild(layout);
        const saved = await grok.dapi.projects.save(project);
        return { projectId: saved.id, layoutId: layout.id };
      });
      projectId = ids.projectId;
      layoutId = ids.layoutId;
      expect(projectId).toBeTruthy();
      expect(layoutId).toBeTruthy();
    },
  );

  await softStep(
    "11. Open the project — TableView and dataframe restore",
    async () => {

      const info = await page.evaluate(async (id) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const proj = await grok.dapi.projects.find(id!);
        await proj.open();
        await new Promise((r) => setTimeout(r, 5000));
        const tv = grok.shell.tv;
        return {
          viewName: grok.shell.v?.name,
          viewType: grok.shell.v?.type,
          rows: tv?.dataFrame?.rowCount,
          viewers: (() => {
            try {
              return Array.from(tv?.viewers ?? []).map((v: any) => v?.type);
            } catch (e) {
              return [];
            }
          })(),
        };
      }, projectId);
      expect(info.viewName).toBe("df_layout");
      expect(info.viewType).toBe("TableView");
      expect(info.rows).toBe(270);
      expect(info.viewers).toContain("Grid");

    },
  );

  await test.step.skip(
    "12. Toolbox > File > Refresh (SKIP — no File pane for script-produced TableViews)",
    async () => {},
  );

  await page.evaluate(
    async (ids) => {
      grok.shell.closeAll();
      try {
        if (ids.projectId) {
          const p = await grok.dapi.projects.find(ids.projectId);
          if (p) await grok.dapi.projects.delete(p);
        }
      } catch (_) {}
      try {
        if (ids.layoutId) {
          const l = await grok.dapi.layouts.find(ids.layoutId);
          if (l) await grok.dapi.layouts.delete(l);
        }
      } catch (_) {}
      const scripts = await grok.dapi.scripts
        .filter('name = "test_Layout"')
        .list();
      for (const sc of scripts)
        try {
          await grok.dapi.scripts.delete(sc);
        } catch (_) {}
    },
    { projectId, layoutId },
  );

  if (stepErrors.length > 0) {
    const summary = stepErrors
      .map((e) => `  - ${e.step}: ${e.error}`)
      .join("\n");
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
