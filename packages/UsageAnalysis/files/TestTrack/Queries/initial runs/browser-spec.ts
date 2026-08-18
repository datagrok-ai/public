import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Queries — browse NorthwindTest and find new_test_query", async ({
  page,
}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  const seedId = await page.evaluate(async () => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    (window as any).grok.shell.windows.showBrowse = true;

    const all = await (window as any).grok.dapi.connections
      .filter("Northwind")
      .list({ pageSize: 30 });
    const conn = all.find(
      (c: any) =>
        c.friendlyName === "NorthwindTest" && c.dataSource === "Postgres",
    );
    let q = await (window as any).grok.dapi.queries
      .filter(`connection.id = "${conn.id}" and name = "new_test_query"`)
      .first()
      .catch(() => null);
    if (!q) {
      const newQ = conn.query("new_test_query", "select * from orders");
      newQ.name = "new_test_query";
      newQ.friendlyName = "new_test_query";
      q = await (window as any).grok.dapi.queries.save(newQ);
    }
    return q.id;
  });

  await page.locator('[name="Browse"]').waitFor({ timeout: 30_000 });

  await softStep("Refresh Browse", async () => {

    await page.goto(`${process.env.DATAGROK_URL}/browse`);
    await page.locator(".d4-tree-view-root").waitFor({ timeout: 15_000 });
    await page.locator('[name="icon-sync"]').click({ timeout: 5_000 });
    await page.locator(".d4-tree-view-root").waitFor({ timeout: 15_000 });
    expect(await page.locator(".d4-tree-view-root").count()).toBeGreaterThan(0);
  });

  await softStep("Browse → Databases → Postgres → NorthwindTest", async () => {

    const dbExpanded = await page.evaluate(async () => {
      const find = (label: string) =>
        Array.from(
          document.querySelectorAll(
            ".d4-tree-view-root .d4-tree-view-group-label",
          ),
        ).find((el) => el.textContent?.trim() === label) as
          | HTMLElement
          | undefined;
      const db = find("Databases");
      if (!db) return { ok: false, stage: "no-db" };
      db.click();
      for (let i = 0; i < 20; i++) {
        if (find("Postgres")) return { ok: true };
        await new Promise((r) => setTimeout(r, 300));
      }
      return { ok: false, stage: "no-pg" };
    });
    expect(dbExpanded.ok).toBe(true);

    const pgExpanded = await page.evaluate(async () => {
      const find = (label: string) =>
        Array.from(
          document.querySelectorAll(
            ".d4-tree-view-root .d4-tree-view-group-label",
          ),
        ).find((el) => el.textContent?.trim() === label) as
          | HTMLElement
          | undefined;
      const pg = find("Postgres");
      if (!pg) return { ok: false, stage: "no-pg" };
      pg.dispatchEvent(
        new MouseEvent("dblclick", { bubbles: true, cancelable: true }),
      );
      for (let i = 0; i < 30; i++) {
        if (find("NorthwindTest")) return { ok: true };
        await new Promise((r) => setTimeout(r, 300));
      }
      return { ok: false, stage: "no-nw" };
    });
    expect(pgExpanded.ok).toBe(true);

    await page.evaluate(async () => {
      const nw = Array.from(
        document.querySelectorAll(
          ".d4-tree-view-root .d4-tree-view-group-label",
        ),
      ).find((el) => el.textContent?.trim() === "NorthwindTest") as HTMLElement;
      nw.dispatchEvent(
        new MouseEvent("dblclick", { bubbles: true, cancelable: true }),
      );
    });
    await page.waitForFunction(
      () => (window as any).grok?.shell?.v?.type === "queries",
      null,
      { timeout: 15_000 },
    );
  });

  await softStep("Type new_test in search field", async () => {

    const typed = await page.evaluate(async () => {
      const input = Array.from(document.querySelectorAll("input")).find(
        (i) =>
          i.offsetParent !== null && i.placeholder?.includes("Search queries"),
      ) as HTMLInputElement | undefined;
      if (!input) return { ok: false, stage: "no-search" };
      input.focus();
      input.value = "new_test";
      input.dispatchEvent(new Event("input", { bubbles: true }));
      input.dispatchEvent(new Event("change", { bubbles: true }));
      input.dispatchEvent(
        new KeyboardEvent("keyup", { bubbles: true, key: "t" }),
      );
      await new Promise((r) => setTimeout(r, 1500));
      const countMatch = document.body.innerText.match(/(\d+)\s*\/\s*(\d+)/);
      return { ok: true, value: input.value, hasCountDisplay: !!countMatch };
    });
    expect(typed.ok).toBe(true);
    expect(typed.value).toBe("new_test");
  });

  await softStep("Context Panel tabs for the query", async () => {
    const tabs = await page.evaluate(async (id) => {
      const q = await (window as any).grok.dapi.queries.find(id);
      (window as any).grok.shell.o = q;
      await new Promise((r) => setTimeout(r, 1500));

      const headers = Array.from(
        document.querySelectorAll(".d4-accordion-pane-header"),
      ).map((h) => (h.textContent ?? "").trim().replace(/\d+$/, ""));
      return [...new Set(headers)];
    }, seedId);

    for (const required of [
      "Details",
      "Run",
      "Query",
      "Transformations",
      "Usage",
      "Sharing",
    ])
      expect(tabs).toContain(required);
  });

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
