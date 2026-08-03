import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Queries — New SQL Query from the products table", async ({ page }) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const w: any = window;
    document.body.classList.add("selenium");
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    w.grok.shell.windows.showBrowse = true;
  });

  // Land on /browse so the standard Databases tree is mounted.
  await page.goto(
    `${process.env.DATAGROK_URL ?? "http://localhost:8888"}/browse`,
  );
  await page.waitForFunction(
    () => {
      return document.querySelectorAll(".d4-tree-view-group-label").length > 5;
    },
    undefined,
    { timeout: 30_000 },
  );

  const queryName = `tt_new_sql_query_${Date.now()}`;

  await softStep(
    "Browse → Databases → Postgres → NorthwindTest → Schemas → public",
    async () => {
      const result = await page.evaluate(async () => {
        const findChildGroupByLabel = (
          parent: ParentNode,
          label: string,
        ): Element | null => {
          const candidates = Array.from(
            parent.querySelectorAll(".d4-tree-view-group-label"),
          ).filter((el) => el.textContent?.trim() === label);
          for (const c of candidates) {
            const g = c.closest(".d4-tree-view-group");
            if (!g) continue;
            if (parent === document || (parent as Element).contains(g))
              return g;
          }
          return null;
        };
        const expandAndWaitFor = async (group: Element, nextLabel: string) => {
          const tri = group.querySelector(
            ":scope > .d4-tree-view-node > .d4-tree-view-tri",
          ) as HTMLElement | null;
          if (tri && !tri.classList.contains("d4-tree-view-tri-expanded"))
            tri.click();
          const deadline = Date.now() + 30_000;
          while (Date.now() < deadline) {
            await new Promise((r) => setTimeout(r, 200));
            const found = findChildGroupByLabel(group, nextLabel);
            if (found) return true;
          }
          return false;
        };
        const path = [
          "Databases",
          "Postgres",
          "NorthwindTest",
          "Schemas",
          "public",
        ];
        let scope: ParentNode = document;
        for (let i = 0; i < path.length; i++) {
          const label = path[i];
          const group = findChildGroupByLabel(scope, label);
          if (!group) return { ok: false, missing: label };
          const next = i + 1 < path.length ? path[i + 1] : "products";
          const ok = await expandAndWaitFor(group, next);
          if (!ok) return { ok: false, expandTimedOutAt: label };
          scope = group;
        }
        return {
          ok: true,
          hasProducts: !!findChildGroupByLabel(scope, "products"),
        };
      });
      expect(result.ok).toBe(true);
      expect(result.hasProducts).toBe(true);
    },
  );

  await softStep(
    "Right-click products → New SQL Query… opens query editor",
    async () => {
      const result = await page.evaluate(async () => {
        const findChildGroupByLabel = (
          parent: ParentNode,
          label: string,
        ): Element | null => {
          const candidates = Array.from(
            parent.querySelectorAll(".d4-tree-view-group-label"),
          ).filter((el) => el.textContent?.trim() === label);
          for (const c of candidates) {
            const g = c.closest(".d4-tree-view-group");
            if (!g) continue;
            if (parent === document || (parent as Element).contains(g))
              return g;
          }
          return null;
        };
        let scope: ParentNode = document;
        for (const label of [
          "Databases",
          "Postgres",
          "NorthwindTest",
          "Schemas",
          "public",
        ]) {
          scope = findChildGroupByLabel(scope, label) as Element;
          if (!scope) return { ok: false, lost: label };
        }
        const productsGroup = findChildGroupByLabel(scope, "products");
        if (!productsGroup) return { ok: false, missing: "products" };
        const tableNode = productsGroup.querySelector(
          ":scope > .d4-tree-view-node",
        ) as HTMLElement | null;
        if (!tableNode) return { ok: false, missing: "tableNode" };
        (
          tableNode.querySelector(
            ".d4-tree-view-group-label",
          ) as HTMLElement | null
        )?.click();
        await new Promise((r) => setTimeout(r, 250));
        tableNode.dispatchEvent(
          new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            button: 2,
          }),
        );
        let target: HTMLElement | undefined;
        const menuDeadline = Date.now() + 8_000;
        while (Date.now() < menuDeadline) {
          await new Promise((r) => setTimeout(r, 150));
          target = Array.from(
            document.querySelectorAll(".d4-menu-popup .d4-menu-item-label"),
          ).find((el) => el.textContent?.trim() === "New SQL Query...") as
            | HTMLElement
            | undefined;
          if (target) break;
        }
        if (!target)
          return { ok: false, missing: "New SQL Query... menu item" };
        target.click();
        // Wait for query editor (CodeMirror + DataQueryView)
        for (let i = 0; i < 80; i++) {
          await new Promise((r) => setTimeout(r, 200));
          const cm: any = document.querySelector(".CodeMirror");
          if (cm && cm.CodeMirror) {
            const v: any = (window as any).grok.shell.v;
            return {
              ok: true,
              viewType: v?.type,
              cm: cm.CodeMirror.getValue(),
            };
          }
        }
        return { ok: false, missing: "CodeMirror editor" };
      });
      expect(result.ok).toBe(true);
      expect(result.viewType).toBe("DataQueryView");
      expect(result.cm).toContain("select * from public.products");
    },
  );

  await softStep(
    "Run via ribbon Play button — inline grid appears",
    async () => {
      const result = await page.evaluate(async () => {
        const playIcon = document.querySelector(
          '[name="icon-play"]',
        ) as HTMLElement | null;
        if (!playIcon) return { ok: false, missing: "icon-play" };
        const ribbonItem = (playIcon.closest(".d4-ribbon-item") ??
          playIcon) as HTMLElement;
        const v: any = (window as any).grok.shell.v;
        const root: HTMLElement | undefined = v?.root;
        if (!root) return { ok: false, missing: "view root" };
        // The Play button click handler is sometimes wired up after a brief delay
        // following editor mount — re-clicking until the inline grid appears.
        const deadline = Date.now() + 60_000;
        let attempts = 0;
        while (Date.now() < deadline) {
          ribbonItem.click();
          playIcon.click();
          attempts++;
          for (let i = 0; i < 12; i++) {
            await new Promise((r) => setTimeout(r, 250));
            const grid = root.querySelector('[name="viewer-Grid"]');
            const canvas = grid?.querySelector("canvas");
            if (canvas) {
              const r2 = grid!.getBoundingClientRect();
              return { ok: true, width: r2.width, height: r2.height, attempts };
            }
          }
        }
        return { ok: false, missing: "inline grid", attempts };
      });
      expect(result.ok).toBe(true);
      expect(result.width).toBeGreaterThan(100);
    },
  );

  await softStep(
    "Run via Toolbox → Actions → Run query… — new view opens",
    async () => {
      const result = await page.evaluate(async () => {
        const w: any = window;
        const before = Array.from(w.grok.shell.views).length;
        const link = Array.from(
          document.querySelectorAll("label.d4-link-action"),
        ).find(
          (el) =>
            el.textContent?.trim() === "Run query..." &&
            (el as HTMLElement).offsetParent !== null,
        ) as HTMLElement | undefined;
        if (!link) return { ok: false, missing: "Run query... link", before };
        link.click();
        for (let i = 0; i < 60; i++) {
          await new Promise((r) => setTimeout(r, 300));
          const after = Array.from(w.grok.shell.views).length;
          if (after > before) {
            // Wait for the result table to settle
            for (let j = 0; j < 30; j++) {
              await new Promise((r) => setTimeout(r, 200));
              const v = w.grok.shell.v;
              if (v?.type === "TableView" && v?.dataFrame?.rowCount > 0)
                return {
                  ok: true,
                  before,
                  after,
                  rowCount: v.dataFrame.rowCount,
                  colCount: v.dataFrame.columns.length,
                  name: v.name,
                };
            }
            return { ok: true, before, after, settled: false };
          }
        }
        return { ok: false, missing: "new view", before };
      });
      expect(result.ok).toBe(true);
      expect(result.after).toBeGreaterThan(result.before);
      expect(result.rowCount).toBeGreaterThan(0);
    },
  );

  await softStep("Save the query", async () => {
    // Switch back to the DataQueryView (the editor)
    await page.evaluate(() => {
      const w: any = window;
      const editor = (Array.from(w.grok.shell.views) as any[]).find(
        (v) => v.type === "DataQueryView",
      );
      if (editor) w.grok.shell.v = editor;
    });
    await page.waitForTimeout(800);
    // Set unique Name in the editor's Name input
    await page.evaluate((name) => {
      const input = Array.from(
        document.querySelectorAll("input.ui-input-editor"),
      ).find(
        (el) =>
          (el as HTMLInputElement).value === "products" &&
          (el as HTMLElement).offsetParent !== null,
      ) as HTMLInputElement | undefined;
      if (!input) throw new Error("no Name input with value=products");
      input.focus();
      document.execCommand("selectAll");
      document.execCommand("insertText", false, name);
    }, queryName);
    await page.waitForTimeout(400);
    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(4000);
    const found = await page.evaluate(async (n) => {
      const q = await (window as any).grok.dapi.queries
        .filter(`name = "${n}"`)
        .first()
        .catch(() => null);
      return !!q;
    }, queryName);
    expect(found).toBe(true);
  });

  // Cleanup — delete the saved query
  await page.evaluate(async (n) => {
    const q = await (window as any).grok.dapi.queries
      .filter(`name = "${n}"`)
      .first()
      .catch(() => null);
    if (q)
      try {
        await (window as any).grok.dapi.queries.delete(q);
      } catch (e) {}
  }, queryName);

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
