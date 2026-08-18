import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

async function expandTreeAndContextMenu(
  page: import("@playwright/test").Page,
  path: string[],
  tableName: string,
  menuLabel: "Get All" | "Get Top 100",
) {

  await page.evaluate(async () => {
    const w: any = window;
    w.grok.shell.windows.showBrowse = true;
    if (document.querySelectorAll(".d4-tree-view-group-label").length < 5) {
      document.querySelector<HTMLElement>('[name="Browse"]')?.click();
      await new Promise((r) => setTimeout(r, 1500));
    }
  });
  await page.waitForFunction(
    () => {
      return document.querySelectorAll(".d4-tree-view-group-label").length > 5;
    },
    undefined,
    { timeout: 30_000 },
  );

  return await page.evaluate(
    async ({ path, tableName, menuLabel }) => {
      const w: any = window;

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
          if (parent === document || (parent as Element).contains(g)) return g;
        }
        return null;
      };

      const expandAndWaitFor = async (
        group: Element,
        nextLabel: string | null,
      ) => {
        const tri = group.querySelector(
          ":scope > .d4-tree-view-node > .d4-tree-view-tri",
        ) as HTMLElement | null;
        if (tri && !tri.classList.contains("d4-tree-view-tri-expanded"))
          tri.click();

        const deadline = Date.now() + 25_000;
        while (Date.now() < deadline) {
          await new Promise((r) => setTimeout(r, 200));
          if (nextLabel) {
            const found = findChildGroupByLabel(group, nextLabel);
            if (found) return true;
          } else {
            const host = group.querySelector(
              ":scope > .d4-tree-view-group-host",
            );
            if (host && host.children.length > 0) return true;
          }
        }
        return false;
      };

      let scope: ParentNode = document;
      for (let i = 0; i < path.length; i++) {
        const label = path[i];
        const group = findChildGroupByLabel(scope, label);
        if (!group) {
          const visibleLabels = Array.from(
            document.querySelectorAll(".d4-tree-view-group-label"),
          )
            .filter((el) => (el as HTMLElement).offsetParent !== null)
            .slice(0, 30)
            .map((el) => el.textContent?.trim());
          return {
            error: `no node "${label}" in path ${path.slice(0, i + 1).join(" > ")}`,
            visibleLabels,
          };
        }
        const nextLabel = i + 1 < path.length ? path[i + 1] : tableName;
        const ok = await expandAndWaitFor(group, nextLabel);
        if (!ok) return { error: `expand timed out at "${label}"` };
        scope = group;
      }

      const tableGroup = findChildGroupByLabel(scope, tableName);
      if (!tableGroup) return { error: `no table "${tableName}"` };
      const tableNode = tableGroup.querySelector(
        ":scope > .d4-tree-view-node",
      ) as HTMLElement | null;
      if (!tableNode) return { error: "no table node" };

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
        ).find((el) => el.textContent?.trim() === menuLabel) as
          | HTMLElement
          | undefined;
        if (target) break;
      }
      if (!target) {
        const items = Array.from(
          document.querySelectorAll(".d4-menu-popup .d4-menu-item-label"),
        ).map((el) => el.textContent?.trim());
        return { error: `no menu item "${menuLabel}"`, items };
      }
      target.click();

      const expectAtLeast = menuLabel === "Get Top 100" ? 100 : 200;
      let last = 0;
      let stable = 0;
      for (let i = 0; i < 120; i++) {
        await new Promise((r) => setTimeout(r, 300));
        const rc = w.grok.shell.tv?.dataFrame?.rowCount ?? 0;
        if (rc >= expectAtLeast && rc === last) {
          stable++;
          if (stable >= 3) break;
        } else {
          last = rc;
          stable = 0;
        }
      }
      const tv = w.grok.shell.tv;
      const df = tv?.dataFrame;
      return {
        viewName: tv?.name,
        rowCount: df?.rowCount,
        colCount: df?.columns?.length,
      };
    },
    { path, tableName, menuLabel },
  );
}

test("Queries — Get All / Get Top 100 (PostgresDart NorthwindTest, Postgres Northwind)", async ({
  page,
}) => {
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
  await page.waitForTimeout(800);

  await softStep(
    "Part 1 — PostgresDart > NorthwindTest > Schemas > public > orders > Get All",
    async () => {
      const r = await expandTreeAndContextMenu(
        page,
        ["Databases", "PostgresDart", "NorthwindTest", "Schemas", "public"],
        "orders",
        "Get All",
      );
      expect((r as any).error, JSON.stringify(r)).toBeUndefined();
      expect((r as any).rowCount).toBeGreaterThan(100);
      expect((r as any).colCount).toBeGreaterThan(0);
    },
  );

  await softStep(
    "Part 1 — PostgresDart > NorthwindTest > Schemas > public > orders > Get Top 100",
    async () => {
      const r = await expandTreeAndContextMenu(
        page,
        ["Databases", "PostgresDart", "NorthwindTest", "Schemas", "public"],
        "orders",
        "Get Top 100",
      );
      expect((r as any).error, JSON.stringify(r)).toBeUndefined();
      expect((r as any).rowCount).toBe(100);
    },
  );

  await softStep(
    "Part 2 — Postgres > Northwind > Schemas > public > orders > Get All",
    async () => {
      const r = await expandTreeAndContextMenu(
        page,
        ["Databases", "Postgres", "Northwind", "Schemas", "public"],
        "orders",
        "Get All",
      );
      expect((r as any).error, JSON.stringify(r)).toBeUndefined();
      expect((r as any).rowCount).toBeGreaterThan(100);
      expect((r as any).colCount).toBeGreaterThan(0);
    },
  );

  await softStep(
    "Part 2 — Postgres > Northwind > Schemas > public > orders > Get Top 100",
    async () => {
      const r = await expandTreeAndContextMenu(
        page,
        ["Databases", "Postgres", "Northwind", "Schemas", "public"],
        "orders",
        "Get Top 100",
      );
      expect((r as any).error, JSON.stringify(r)).toBeUndefined();
      expect((r as any).rowCount).toBe(100);
    },
  );

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
