import { test, expect } from "@playwright/test";
import {
  baseUrl,
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

const CONN_NAME = "PostgreSQLDBTests2";
const SERVER = "db.datagrok.ai";
const PORT = "54327";
const DB = "test";
const LOGIN = "datagrok";
const PASSWORD =
  process.env.DATAGROK_PG_TEST_PASSWORD ??
  "WZNFYTSDwu8TTfN6RQ2Yp5VKvbJD0Fddslpm";

const QUERIES: { name: string; sql: string }[] = [
  {
    name: "TestCreateTable",
    sql: "CREATE TABLE tmp_table_test (id bigint, name varchar)",
  },
  {
    name: "TestInsertData",
    sql: "INSERT INTO tmp_table_test VALUES (1, 'test')",
  },
  {
    name: "TestUpdateData",
    sql: "UPDATE tmp_table_test SET name = 'test' WHERE id = 1",
  },
  { name: "TestDropTable", sql: "DROP TABLE tmp_table_test" },
];

test("Connections / External Provider", async ({ page }) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const existing = await (window as any).grok.dapi.connections
      .filter(`name = "PostgreSQLDBTests2"`)
      .list();
    for (const c of existing) {
      try {
        const qs = await (window as any).grok.dapi.queries
          .filter(`connection.id = "${c.id}"`)
          .list();
        for (const q of qs)
          try {
            await (window as any).grok.dapi.queries.delete(q);
          } catch (e) {}
        await (window as any).grok.dapi.connections.delete(c);
      } catch (e) {}
    }
    for (const n of [
      "TestCreateTable",
      "TestInsertData",
      "TestUpdateData",
      "TestDropTable",
    ]) {
      const list = await (window as any).grok.dapi.queries
        .filter(`name = "${n}"`)
        .list();
      for (const q of list)
        try {
          await (window as any).grok.dapi.queries.delete(q);
        } catch (e) {}
    }
  });

  await page.goto(`${baseUrl}/browse`, { timeout: 30_000 });
  await page.waitForTimeout(5000);

  await softStep("Step 1: Browse > Databases > Postgres", async () => {
    await page.evaluate(async () => {
      const findLabel = (t: string) =>
        Array.from(
          document.querySelectorAll(
            ".d4-tree-view-group-label, .d4-tree-view-item-label",
          ),
        ).find((e) => e.textContent?.trim() === t) as HTMLElement | undefined;
      let dbs = findLabel("Databases");
      dbs?.click();
      await new Promise((r) => setTimeout(r, 800));
      let pg = findLabel("Postgres");
      if (!pg) {
        dbs?.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
        await new Promise((r) => setTimeout(r, 1500));
        pg = findLabel("Postgres");
      }
      pg?.click();
      await new Promise((r) => setTimeout(r, 800));
    });
    await page.waitForFunction(
      () => {
        return Array.from(
          document.querySelectorAll(
            ".d4-tree-view-group-label, .d4-tree-view-item-label",
          ),
        ).some((el) => el.textContent?.trim() === "Postgres");
      },
      null,
      { timeout: 30_000 },
    );
    const found = await page.evaluate(() => {
      return Array.from(
        document.querySelectorAll(
          ".d4-tree-view-group-label, .d4-tree-view-item-label",
        ),
      ).some((el) => el.textContent?.trim() === "Postgres");
    });
    expect(found).toBe(true);
  });

  await softStep("Step 2: Open Add new connection dialog", async () => {
    await page.evaluate(async () => {
      const labels = Array.from(
        document.querySelectorAll(
          ".d4-tree-view-group-label, .d4-tree-view-item-label",
        ),
      ).filter((el) => el.textContent?.trim() === "Postgres");
      const node =
        labels
          .map(
            (el) =>
              el.closest(
                ".d4-tree-view-node, .d4-tree-view-group",
              ) as HTMLElement | null,
          )
          .find((n) => !!n && (n as HTMLElement).offsetParent !== null) ??
        labels[0]?.closest(".d4-tree-view-node, .d4-tree-view-group");
      (node as HTMLElement | null)?.dispatchEvent(
        new MouseEvent("contextmenu", {
          bubbles: true,
          cancelable: true,
          button: 2,
        }),
      );
      await new Promise((r) => setTimeout(r, 600));
      const item = Array.from(
        document.querySelectorAll(".d4-menu-popup .d4-menu-item-label"),
      ).find((el) => el.textContent?.trim() === "New connection...");
      (item as HTMLElement | undefined)?.click();
    });
    await page.waitForFunction(
      () => !!document.querySelector(".d4-dialog"),
      null,
      { timeout: 15_000 },
    );
    const hasDialog = await page.evaluate(
      () => !!document.querySelector(".d4-dialog"),
    );
    expect(hasDialog).toBe(true);
  });

  await softStep(
    "Step 3: Fill connection form (Name/Server/Port/Db/Login/Password) and Save",
    async () => {
      await page.waitForFunction(
        () => !!document.querySelector('input[name="input-Name"]'),
        null,
        { timeout: 60_000 },
      );
      await page.waitForFunction(
        () => !!document.querySelector('input[name="input-Password"]'),
        null,
        { timeout: 30_000 },
      );
      await page.waitForTimeout(500);
      await page.evaluate(
        ({ name, server, port, db, login, pwd }) => {
          const setText = (selector: string, value: string) => {
            const input = Array.from(document.querySelectorAll(selector)).find(
              (el) => (el as HTMLElement).offsetParent !== null,
            ) as HTMLInputElement | undefined;
            if (!input) return false;
            const setter = Object.getOwnPropertyDescriptor(
              window.HTMLInputElement.prototype,
              "value",
            )!.set!;
            input.focus();
            setter.call(input, "");
            input.dispatchEvent(new Event("input", { bubbles: true }));
            setter.call(input, value);
            input.dispatchEvent(new Event("input", { bubbles: true }));
            input.dispatchEvent(new Event("change", { bubbles: true }));
            input.dispatchEvent(
              new KeyboardEvent("keydown", { bubbles: true, key: "Tab" }),
            );
            input.blur();
            return true;
          };
          setText('input[name="input-Name"]', name);
          setText('input[name="input-Server"]', server);
          setText('input[name="input-Port"]', port);
          setText('input[name="input-Db"]', db);
          setText('input[name="input-Login"]', login);
          setText('input[name="input-Password"]', pwd);
        },
        {
          name: CONN_NAME,
          server: SERVER,
          port: PORT,
          db: DB,
          login: LOGIN,
          pwd: PASSWORD,
        },
      );
      await page.waitForTimeout(500);

      await page.evaluate(() => {
        const dialog = document.querySelector(".d4-dialog")!;
        const ok =
          (dialog.querySelector('[name="button-OK"]') as HTMLElement) ??
          (Array.from(dialog.querySelectorAll("button")).find(
            (b) => b.textContent?.trim() === "OK",
          ) as HTMLElement);
        ok?.click();
      });
      await page.waitForTimeout(2500);

      await page.waitForFunction(
        async () => {
          const c = await (window as any).grok.dapi.connections
            .filter('name = "PostgreSQLDBTests2"')
            .first();
          return !!c;
        },
        null,
        { timeout: 60_000 },
      );

      const exists = await page.evaluate(async () => {
        const c = await (window as any).grok.dapi.connections
          .filter('name = "PostgreSQLDBTests2"')
          .first();
        return !!c;
      });
      expect(exists).toBe(true);
    },
  );

  for (const { name: qName, sql } of QUERIES) {
    await softStep(`Step 4.${qName}: create + run query`, async () => {
      await page.goto(`${baseUrl}/browse`, { timeout: 30_000 });
      await page.waitForTimeout(3000);
      await page.evaluate(async () => {
        const findLabel = (t: string) =>
          Array.from(
            document.querySelectorAll(
              ".d4-tree-view-group-label, .d4-tree-view-item-label",
            ),
          ).find((e) => e.textContent?.trim() === t) as HTMLElement | undefined;
        let dbs = findLabel("Databases");
        dbs?.click();
        await new Promise((r) => setTimeout(r, 800));
        let pg = findLabel("Postgres");
        if (!pg) {
          dbs?.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
          await new Promise((r) => setTimeout(r, 1500));
          pg = findLabel("Postgres");
        }
        pg?.click();
        await new Promise((r) => setTimeout(r, 1500));
      });
      await page.waitForFunction(
        (connName) => {
          return Array.from(
            document.querySelectorAll(
              ".d4-tree-view-group-label, .d4-tree-view-item-label",
            ),
          ).some((el) => el.textContent?.trim() === connName);
        },
        CONN_NAME,
        { timeout: 30_000 },
      );

      await page.evaluate(async (connName) => {
        const labels = Array.from(
          document.querySelectorAll(
            ".d4-tree-view-group-label, .d4-tree-view-item-label",
          ),
        ).filter((el) => el.textContent?.trim() === connName);
        const node =
          labels
            .map(
              (el) =>
                el.closest(
                  ".d4-tree-view-node, .d4-tree-view-group",
                ) as HTMLElement | null,
            )
            .find((n) => !!n && (n as HTMLElement).offsetParent !== null) ??
          labels[0]?.closest(".d4-tree-view-node, .d4-tree-view-group");
        (node as HTMLElement | null)?.dispatchEvent(
          new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            button: 2,
          }),
        );
        await new Promise((r) => setTimeout(r, 600));
        const item = Array.from(
          document.querySelectorAll(".d4-menu-popup .d4-menu-item-label"),
        ).find((el) => el.textContent?.trim() === "New Query...");
        (item as HTMLElement | undefined)?.click();
      }, CONN_NAME);

      await page.waitForFunction(
        () => !!document.querySelector(".CodeMirror"),
        null,
        { timeout: 30_000 },
      );
      await page.waitForTimeout(700);

      await page.evaluate(
        ({ queryName, body }) => {
          const setter = Object.getOwnPropertyDescriptor(
            window.HTMLInputElement.prototype,
            "value",
          )!.set!;
          const nameInput = Array.from(
            document.querySelectorAll('input[name="input-Name"]'),
          ).find((el) => (el as HTMLElement).offsetParent !== null) as
            | HTMLInputElement
            | undefined;
          if (nameInput) {
            nameInput.focus();
            setter.call(nameInput, "");
            nameInput.dispatchEvent(new Event("input", { bubbles: true }));
            setter.call(nameInput, queryName);
            nameInput.dispatchEvent(new Event("input", { bubbles: true }));
            nameInput.dispatchEvent(new Event("change", { bubbles: true }));
            nameInput.blur();
          }
          const cm = document.querySelector(".CodeMirror") as any;
          cm?.CodeMirror?.setValue(body);
        },
        { queryName: qName, body: sql },
      );

      await page.waitForTimeout(400);
      await page.evaluate(() =>
        (
          document.querySelector('[name="button-Save"]') as HTMLElement | null
        )?.click(),
      );

      await page.waitForFunction(
        async (queryName) => {
          const q = await (window as any).grok.dapi.queries
            .filter(`name = "${queryName}"`)
            .first();
          return !!q;
        },
        qName,
        { timeout: 30_000 },
      );

      const ran = await page.evaluate(async (queryName) => {
        const q = await (window as any).grok.dapi.queries
          .filter(`name = "${queryName}"`)
          .first();
        try {
          await q.executeTable();
          return true;
        } catch (e) {
          return false;
        }
      }, qName);
      expect(ran).toBe(true);
    });
  }

  await softStep("Step 5: Delete PostgreSQLDBTests2 connection", async () => {
    await page.goto(`${baseUrl}/browse`, { timeout: 30_000 });
    await page.waitForTimeout(3000);
    await page.evaluate(async () => {
      const findLabel = (t: string) =>
        Array.from(
          document.querySelectorAll(
            ".d4-tree-view-group-label, .d4-tree-view-item-label",
          ),
        ).find((e) => e.textContent?.trim() === t) as HTMLElement | undefined;
      let dbs = findLabel("Databases");
      dbs?.click();
      await new Promise((r) => setTimeout(r, 800));
      let pg = findLabel("Postgres");
      if (!pg) {
        dbs?.dispatchEvent(new MouseEvent("dblclick", { bubbles: true }));
        await new Promise((r) => setTimeout(r, 1500));
        pg = findLabel("Postgres");
      }
      pg?.click();
      await new Promise((r) => setTimeout(r, 1500));
    });
    await page.waitForFunction(
      (connName) => {
        return Array.from(
          document.querySelectorAll(
            ".d4-tree-view-group-label, .d4-tree-view-item-label",
          ),
        ).some((el) => el.textContent?.trim() === connName);
      },
      CONN_NAME,
      { timeout: 30_000 },
    );

    await page.evaluate(async (connName) => {
      const node = Array.from(
        document.querySelectorAll(
          ".d4-tree-view-group-label, .d4-tree-view-item-label",
        ),
      ).find((el) => el.textContent?.trim() === connName);
      (
        node?.closest(".d4-tree-view-node, .d4-tree-view-group") as HTMLElement
      )?.dispatchEvent(
        new MouseEvent("contextmenu", {
          bubbles: true,
          cancelable: true,
          button: 2,
        }),
      );
      await new Promise((r) => setTimeout(r, 500));
      const item = Array.from(
        document.querySelectorAll(".d4-menu-popup .d4-menu-item-label"),
      ).find((el) => el.textContent?.trim() === "Delete...");
      (item as HTMLElement | undefined)?.click();
    }, CONN_NAME);

    await page.waitForFunction(
      () => !!document.querySelector(".d4-dialog"),
      null,
      { timeout: 10_000 },
    );
    await page.evaluate(() => {
      const dialog = document.querySelector(".d4-dialog")!;
      const del =
        (dialog.querySelector('[name="button-DELETE"]') as HTMLElement) ??
        (Array.from(dialog.querySelectorAll("button")).find(
          (b) => b.textContent?.trim() === "DELETE",
        ) as HTMLElement);
      del?.click();
    });

    await page.waitForFunction(
      async () => {
        const c = await (window as any).grok.dapi.connections
          .filter('name = "PostgreSQLDBTests2"')
          .first();
        return !c;
      },
      null,
      { timeout: 30_000 },
    );

    const gone = await page.evaluate(async () => {
      const c = await (window as any).grok.dapi.connections
        .filter('name = "PostgreSQLDBTests2"')
        .first();
      return !c;
    });
    expect(gone).toBe(true);
  });

  await page.evaluate(async () => {
    for (const n of [
      "TestCreateTable",
      "TestInsertData",
      "TestUpdateData",
      "TestDropTable",
    ]) {
      const list = await (window as any).grok.dapi.queries
        .filter(`name = "${n}"`)
        .list();
      for (const q of list)
        try {
          await (window as any).grok.dapi.queries.delete(q);
        } catch (e) {}
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors
      .map((e) => `  - ${e.step}: ${e.error}`)
      .join("\n");
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
