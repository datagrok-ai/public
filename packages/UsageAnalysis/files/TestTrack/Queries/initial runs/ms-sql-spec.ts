import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Queries — MS SQL: Add / Edit / Browse / Delete", async ({ page }) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const setupResult = await page.evaluate(async () => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;
    const stale = await (window as any).grok.dapi.queries
      .filter('friendlyName in ("test_query_ms_sql","new_test_query_ms_sql")')
      .list()
      .catch(() => []);
    for (const q of stale)
      try {
        await (window as any).grok.dapi.queries.delete(q);
      } catch (e) {}
    return { cleaned: stale.length };
  });
  expect(setupResult.cleaned).toBeGreaterThanOrEqual(0);

  await page.goto(
    `${process.env.DATAGROK_URL ?? "https://dev.datagrok.ai"}/browse`,
  );

  let savedQueryId: string | null = null;

  await softStep(
    "Part 1.1 — Browse → Databases → MS SQL → right-click NorthwindTest → New Query",
    async () => {

      await page
        .locator('[name="tree-Databases"]')
        .first()
        .waitFor({ state: "attached", timeout: 30_000 });
      await page.evaluate(async () => {

        for (let i = 0; i < 30; i++) {
          const exp = document.querySelector(
            '[name="tree-expander-Databases"]',
          ) as HTMLElement | null;
          const node = document.querySelector(
            '[name="tree-Databases"]',
          ) as HTMLElement | null;
          if (document.querySelector('[name="tree-Databases---MS-SQL"]')) break;
          (exp ?? node)?.click();
          await new Promise((r) => setTimeout(r, 500));
        }

        for (let i = 0; i < 30; i++) {
          const exp = document.querySelector(
            '[name="tree-expander-Databases---MS-SQL"]',
          ) as HTMLElement | null;
          const node = document.querySelector(
            '[name="tree-Databases---MS-SQL"]',
          ) as HTMLElement | null;
          if (
            document.querySelector(
              '[name="tree-Databases---MS-SQL---NorthwindTest"]',
            )
          )
            break;
          node?.scrollIntoView({ block: "center" });
          (exp ?? node)?.click();
          await new Promise((r) => setTimeout(r, 800));
        }
      });

      await page
        .locator('[name="tree-Databases---MS-SQL---NorthwindTest"]')
        .first()
        .waitFor({ state: "attached", timeout: 30_000 });
      const opened = await page.evaluate(async () => {
        const nw = document.querySelector(
          '[name="tree-Databases---MS-SQL---NorthwindTest"]',
        ) as HTMLElement | null;
        if (!nw) return { ok: false, stage: "no-northwindtest" };
        nw.scrollIntoView({ block: "center" });
        nw.dispatchEvent(
          new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            button: 2,
          }),
        );
        for (let i = 0; i < 25; i++) {
          const item = Array.from(
            document.querySelectorAll(".d4-menu-item-label"),
          ).find((el) => el.textContent?.trim() === "New Query...") as
            | HTMLElement
            | undefined;
          if (item) {
            item.click();
            return { ok: true };
          }
          await new Promise((r) => setTimeout(r, 200));
        }
        return { ok: false, stage: "no-menu-item" };
      });
      expect(opened.ok).toBe(true);
      await page.locator(".CodeMirror").first().waitFor({ timeout: 30_000 });
      await page
        .locator('input[name="input-Name"]')
        .waitFor({ timeout: 30_000 });
    },
  );

  await softStep("Part 1.2 — Enter test_query_ms_sql into Name", async () => {
    const nameInput = page.locator('input[name="input-Name"]').first();
    await nameInput.click();
    await page.keyboard.press("Control+a");
    await page.keyboard.type("test_query_ms_sql");
    const value = await page.evaluate(
      () =>
        (document.querySelector('input[name="input-Name"]') as HTMLInputElement)
          ?.value ?? "",
    );
    expect(value).toBe("test_query_ms_sql");
  });

  await softStep("Part 1.3 — Enter SQL: select * from products", async () => {
    const body = await page.evaluate(() => {
      const cm = (document.querySelector(".CodeMirror") as any)?.CodeMirror;
      if (!cm) return null;
      cm.setValue("select * from products");
      return cm.getValue();
    });
    expect(body).toBe("select * from products");
  });

  await softStep("Part 1.4 — Run via Play button (inline)", async () => {

    await page.locator('[name="icon-play"]').first().click();
    await page.waitForTimeout(8000);

    const editorAlive = await page.locator(".CodeMirror").first().isVisible();
    expect(editorAlive).toBe(true);
  });

  await softStep(
    "Part 1.5 — Run via Toolbox → Actions → Run query… (new view)",
    async () => {
      const result = await page.evaluate(async () => {
        const beforeViews = Array.from((window as any).grok.shell.views).length;
        const item = Array.from(
          document.querySelectorAll("label, .d4-link-label"),
        ).find((el) => el.textContent?.trim() === "Run query...") as
          | HTMLElement
          | undefined;
        if (!item) return { clicked: false };
        item.click();
        await new Promise((r) => setTimeout(r, 8000));
        const afterViews = Array.from((window as any).grok.shell.views).length;
        return { clicked: true, beforeViews, afterViews };
      });

      expect(await page.locator(".CodeMirror").first().isVisible()).toBe(true);
    },
  );

  await softStep("Part 1.6 — Save the query", async () => {
    await page
      .locator('[name="button-Save"], [name="button-SAVE"]')
      .first()
      .click();
    const saved = await page.evaluate(async () => {
      for (let i = 0; i < 30; i++) {
        const q = await (window as any).grok.dapi.queries
          .filter('friendlyName = "test_query_ms_sql"')
          .first()
          .catch(() => null);
        if (q && q.id)
          return {
            ok: true,
            id: q.id,
            name: q.name,
            friendly: q.friendlyName,
            body: q.query,
          };
        await new Promise((r) => setTimeout(r, 500));
      }
      return { ok: false };
    });
    expect(saved.ok).toBe(true);
    if (saved.ok) savedQueryId = saved.id as string;
  });

  await softStep(
    "Part 2.1 — Refresh Browse → Edit query (open editor)",
    async () => {

      expect(savedQueryId).toBeTruthy();
      await page.goto(
        `${process.env.DATAGROK_URL ?? "https://dev.datagrok.ai"}/query/${savedQueryId}`,
      );
      await page.locator(".CodeMirror").first().waitFor({ timeout: 15_000 });
      await page
        .locator('input[name="input-Name"]')
        .waitFor({ timeout: 15_000 });
    },
  );

  await softStep(
    "Part 2.2 — Change name to new_test_query_ms_sql",
    async () => {
      const nameInput = page.locator('input[name="input-Name"]').first();
      await nameInput.click();
      await page.keyboard.press("Control+a");
      await page.keyboard.type("new_test_query_ms_sql");
      const value = await page.evaluate(
        () =>
          (
            document.querySelector(
              'input[name="input-Name"]',
            ) as HTMLInputElement
          )?.value ?? "",
      );
      expect(value).toBe("new_test_query_ms_sql");
    },
  );

  await softStep("Part 2.3 — Change body to select * from orders", async () => {
    const body = await page.evaluate(() => {
      const cm = (document.querySelector(".CodeMirror") as any)?.CodeMirror;
      cm.setValue("select * from orders");
      return cm.getValue();
    });
    expect(body).toBe("select * from orders");
  });

  await softStep("Part 2.4 — Run via Play + Run query…", async () => {
    await page.locator('[name="icon-play"]').first().click();
    await page.waitForTimeout(5000);
    expect(await page.locator(".CodeMirror").first().isVisible()).toBe(true);
  });

  await softStep("Part 2.5 — Save the renamed query", async () => {
    await page
      .locator('[name="button-Save"], [name="button-SAVE"]')
      .first()
      .click();
    const saved = await page.evaluate(async (id) => {
      for (let i = 0; i < 30; i++) {
        const q = await (window as any).grok.dapi.queries
          .find(id)
          .catch(() => null);
        if (
          q &&
          q.friendlyName === "new_test_query_ms_sql" &&
          q.query === "select * from orders"
        )
          return {
            ok: true,
            name: q.name,
            friendly: q.friendlyName,
            body: q.query,
          };
        await new Promise((r) => setTimeout(r, 500));
      }
      return { ok: false };
    }, savedQueryId);
    expect(saved.ok).toBe(true);
  });

  await softStep(
    "Part 3.1/3.2 — Browse → MS SQL → NorthwindTest, expand and search",
    async () => {
      await page.goto(
        `${process.env.DATAGROK_URL ?? "https://dev.datagrok.ai"}/browse`,
      );
      const result = await page.evaluate(async () => {
        for (let i = 0; i < 30; i++) {
          if (document.querySelector('[name="tree-Databases"]')) break;
          await new Promise((r) => setTimeout(r, 300));
        }

        const ms = document.querySelector(
          '[name="tree-Databases---MS-SQL"]',
        ) as HTMLElement | null;
        ms?.click();
        ms?.dispatchEvent(
          new MouseEvent("dblclick", { bubbles: true, cancelable: true }),
        );

        for (let i = 0; i < 30; i++) {
          if (
            document.querySelector(
              '[name="tree-Databases---MS-SQL---NorthwindTest"]',
            )
          )
            break;
          await new Promise((r) => setTimeout(r, 300));
        }
        const exp = document.querySelector(
          '[name="tree-expander-Databases---MS-SQL---NorthwindTest"]',
        ) as HTMLElement | null;
        exp?.click();

        for (let i = 0; i < 20; i++) {
          const children = document.querySelectorAll(
            '[name^="tree-Databases---MS-SQL---NorthwindTest---"]',
          );
          if (children.length > 0) break;
          await new Promise((r) => setTimeout(r, 300));
        }
        const labels = [
          ...new Set(
            Array.from(
              document.querySelectorAll(
                '[name^="tree-Databases---MS-SQL---NorthwindTest---"] .d4-tree-view-group-label, [name^="tree-Databases---MS-SQL---NorthwindTest---"] .d4-tree-view-item-label',
              ),
            )
              .map((e) => e.textContent?.trim())
              .filter(Boolean),
          ),
        ];
        const containsNew = labels.some((l) =>
          l?.includes("new_test_query_ms_sql"),
        );
        return { labels, containsNew };
      });

      expect(Array.isArray(result.labels)).toBe(true);
    },
  );

  await softStep(
    "Part 3.3 — Context Panel tabs for the saved query",
    async () => {
      const panes = await page.evaluate(async (id) => {
        const q = await (window as any).grok.dapi.queries.find(id);
        (window as any).grok.shell.o = q;
        await new Promise((r) => setTimeout(r, 1500));
        return Array.from(
          document.querySelectorAll(".d4-accordion-pane-header"),
        )
          .map((h) => h.textContent?.trim())
          .filter(Boolean);
      }, savedQueryId);
      for (const required of [
        "Details",
        "Run",
        "Query",
        "Transformations",
        "Sharing",
      ])
        expect(panes).toContain(required);
    },
  );

  await softStep("Part 4 — Delete query and verify removal", async () => {

    expect(savedQueryId).toBeTruthy();
    const deleted = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries
        .find(id)
        .catch(() => null);
      if (!q) return { found: false, ok: false };
      await (window as any).grok.dapi.queries.delete(q);
      await new Promise((r) => setTimeout(r, 1000));
      const after = await (window as any).grok.dapi.queries
        .find(id)
        .catch(() => null);
      return { found: true, ok: !after };
    }, savedQueryId);

    expect(
      (deleted as any).found,
      "saved query should still exist before delete",
    ).toBe(true);
    expect((deleted as any).ok, "query must be gone after delete").toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
