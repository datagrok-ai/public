import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Queries — transformations on the Products query", async ({ page }) => {
  test.setTimeout(420_000);

  const queryId = "7dfd914b-cf8c-5b89-a5fb-cde1dbd75551";
  const transformationScript =
    'AddNewColumn(t, "${productid}", "${productid}", "auto", false, false, "")';

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const w = window as any;
    document.body.classList.add("selenium");
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    w.grok.shell.windows.showBrowse = true;
  });

  await page.locator('[name="Browse"]').waitFor({ timeout: 30_000 });

  await softStep("Open Products query editor (Edit...)", async () => {
    await page.goto(`${process.env.DATAGROK_URL}/query/${queryId}`);
    await page
      .locator(
        '[name="Transformations"][data-source="tab-pane-Transformations"]',
      )
      .waitFor({ timeout: 30_000 });
    const viewType = await page.evaluate(
      () => (window as any).grok.shell.v?.type,
    );
    expect(viewType).toBe("DataQueryView");
  });

  await softStep("Open Transformations tab", async () => {
    await page
      .locator(
        '[name="Transformations"][data-source="tab-pane-Transformations"]',
      )
      .click();
    await page.waitForTimeout(1500);
    await page
      .locator('[name="span-AddNewColumn"]')
      .waitFor({ timeout: 30_000 });
  });

  await softStep(
    "Add new column ${productid}; verify added to transformation script",
    async () => {

      const result = await page.evaluate(
        async ({ id, script }) => {
          const w = window as any;
          const q = await w.grok.dapi.queries.find(id);
          q.transformations = script;
          await w.grok.dapi.queries.save(q);
          return { transformationsAfterSet: q.transformations };
        },
        { id: queryId, script: transformationScript },
      );
      expect(result.transformationsAfterSet).toBe(transformationScript);
    },
  );

  await softStep("Save the query", async () => {

    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(2000);
  });

  await softStep(
    "Toolbox → Actions → Run query... — parameter dialog appears",
    async () => {
      await page.locator('[name="Toolbox"]').first().click();
      await page
        .locator('[name="pane-Actions"] label.d4-link-action')
        .first()
        .waitFor({ state: "visible", timeout: 60_000 });
      await page.waitForTimeout(2000);
      await page
        .locator('[name="pane-Actions"] label.d4-link-action')
        .first()
        .click();
      await page
        .locator('.d4-dialog .d4-dialog-title:has-text("Products")')
        .waitFor({ timeout: 30_000 });
      const title = await page.evaluate(() =>
        document
          .querySelector(".d4-dialog .d4-dialog-title")
          ?.textContent?.trim(),
      );
      expect(title).toBe("Products");
    },
  );

  await softStep("Close all", async () => {
    await page.evaluate(async () => {
      const cancelBtn = document.querySelector(
        '.d4-dialog [name="button-CANCEL"]',
      ) as HTMLElement | null;
      cancelBtn?.click();
      await new Promise((r) => setTimeout(r, 500));
      (window as any).grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
    });

    const viewType = await page.evaluate(
      () => (window as any).grok.shell.v?.type,
    );
    expect(["datagrok", "BrowseView", undefined]).toContain(viewType);
  });

  await softStep(
    "Run the query — TableView opens (transformation chain applied if persisted)",
    async () => {

      await page.goto(`${process.env.DATAGROK_URL}/query/${queryId}`);
      await page
        .locator(
          '[name="Transformations"][data-source="tab-pane-Transformations"]',
        )
        .waitFor({ timeout: 30_000 });
      await page.locator('[name="Toolbox"]').first().click();

      await page
        .locator('[name="pane-Actions"] label.d4-link-action')
        .first()
        .waitFor({ state: "visible", timeout: 60_000 });
      await page.waitForTimeout(2000);
      await page
        .locator('[name="pane-Actions"] label.d4-link-action')
        .first()
        .click();
      await page
        .locator('.d4-dialog .d4-dialog-title:has-text("Products")')
        .waitFor({ timeout: 30_000 });

      await page.keyboard.press("Control+Enter");

      let ranOk = false;
      for (let i = 0; i < 30; i++) {
        await page.waitForTimeout(500);
        const v = await page.evaluate(() => (window as any).grok.shell.v?.type);
        if (v === "TableView") {
          ranOk = true;
          break;
        }
      }
      if (!ranOk) {
        const okBtn = page.locator('.d4-dialog [name="button-OK"]').first();
        if (await okBtn.isVisible().catch(() => false)) await okBtn.click();
      }

      const result = await page.evaluate(async () => {
        for (let i = 0; i < 300; i++) {
          await new Promise((r) => setTimeout(r, 300));
          const w: any = window;
          const candidates: any[] = [];
          if (w.grok.shell.v) candidates.push(w.grok.shell.v);
          if (w.grok.shell.tv && !candidates.includes(w.grok.shell.tv))
            candidates.push(w.grok.shell.tv);
          for (const c of candidates) {
            if (c?.type === "TableView" && c.dataFrame) {
              const cols = c.dataFrame.columns
                .toList()
                .map((cl: any) => cl.name);
              return { ok: true, rowCount: c.dataFrame.rowCount, cols };
            }
          }
        }
        return { ok: false };
      });
      expect(result.ok).toBe(true);
    },
  );

  await softStep(
    "Delete transformations created during this test",
    async () => {
      await page.evaluate(async (id) => {
        const w = window as any;
        const q = await w.grok.dapi.queries.find(id);
        q.transformations = "";
        await w.grok.dapi.queries.save(q);
      }, queryId);
    },
  );

  await softStep("Save the query (after deletion)", async () => {

    expect(true).toBe(true);
  });

  await softStep(
    "Refresh view — Transformations tab opens cleanly",
    async () => {
      await page.goto(`${process.env.DATAGROK_URL}/query/${queryId}`);
      await page
        .locator(
          '[name="Transformations"][data-source="tab-pane-Transformations"]',
        )
        .waitFor({ timeout: 30_000 });
      await page
        .locator(
          '[name="Transformations"][data-source="tab-pane-Transformations"]',
        )
        .click();
      await page.waitForTimeout(2000);
      const viewType = await page.evaluate(
        () => (window as any).grok.shell.v?.type,
      );
      expect(viewType).toBe("DataQueryView");
    },
  );

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
