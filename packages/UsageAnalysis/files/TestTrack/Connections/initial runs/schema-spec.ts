import { test, expect } from "@playwright/test";
import {
  baseUrl,
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Connections / Schema", async ({ page }) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await page.goto(`${baseUrl}/connect?browse=connections`, {
    waitUntil: "networkidle",
    timeout: 30000,
  });
  await page.evaluate(() => {
    document.body.classList.add("selenium");
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });

  await softStep(
    "Step 1-2: Browse > Databases > Postgres expanded, Northwind visible",
    async () => {
      await page.evaluate(async () => {
        for (let i = 0; i < 60; i++) {
          const found = Array.from(
            document.querySelectorAll(".d4-tree-view-group-label"),
          ).some((el) => el.textContent?.trim() === "Postgres");
          if (found) break;
          await new Promise((r) => setTimeout(r, 500));
        }
        const postgres = Array.from(
          document.querySelectorAll(".d4-tree-view-group-label"),
        ).find((el) => el.textContent?.trim() === "Postgres") as
          | HTMLElement
          | undefined;
        const row = postgres?.closest(".d4-tree-view-node");
        const tri = row?.querySelector(
          ".d4-tree-view-tri, .d4-tree-view-group-tri",
        ) as HTMLElement | null;
        tri?.click();
        for (let i = 0; i < 30; i++) {
          const visible = Array.from(
            document.querySelectorAll(".d4-tree-view-group-label"),
          ).some((el) => el.textContent?.trim() === "Northwind");
          if (visible) break;
          await new Promise((r) => setTimeout(r, 300));
        }
      });
      const northwindVisible = await page.evaluate(() =>
        Array.from(document.querySelectorAll(".d4-tree-view-group-label")).some(
          (el) => el.textContent?.trim() === "Northwind",
        ),
      );
      expect(northwindVisible).toBe(true);
    },
  );

  await softStep(
    "Step 3: Right-click Northwind > Browse opens connection view",
    async () => {
      const result = await page.evaluate(async () => {
        const postgres = Array.from(
          document.querySelectorAll(".d4-tree-view-group-label"),
        ).find((el) => el.textContent?.trim() === "Postgres");
        const postgresHost = postgres
          ?.closest(".d4-tree-view-group")
          ?.querySelector(":scope > .d4-tree-view-group-host");
        const groups = postgresHost
          ? Array.from(postgresHost.children).filter((c) =>
              c.classList.contains("d4-tree-view-group"),
            )
          : [];
        const northwindGroup = groups.find(
          (g) =>
            g
              .querySelector(
                ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
              )
              ?.textContent?.trim() === "Northwind",
        );
        const nRow = northwindGroup?.querySelector(
          ":scope > .d4-tree-view-node",
        ) as HTMLElement | null;
        nRow?.scrollIntoView({ block: "center" });
        nRow?.dispatchEvent(
          new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            button: 2,
          }),
        );
        await new Promise((r) => setTimeout(r, 700));
        const browseItem = document.querySelector(
          '[name="div-Browse"]',
        ) as HTMLElement | null;
        (browseItem as HTMLElement | null)?.click();
        await new Promise((r) => setTimeout(r, 3500));
        return {
          browseClicked: !!browseItem,
          viewName: (window as any).grok?.shell?.v?.name,
          viewType: (window as any).grok?.shell?.v?.type,
        };
      });
      expect(result.browseClicked).toBe(true);
      expect(result.viewName).toBe("PostgresNorthwind");
    },
  );

  await softStep(
    "Step 4: Tree exposes Schemas > public > customers with DB-table menu",
    async () => {
      const result = await page.evaluate(async () => {
        const postgres = Array.from(
          document.querySelectorAll(".d4-tree-view-group-label"),
        ).find((el) => el.textContent?.trim() === "Postgres");
        const postgresHost = postgres
          ?.closest(".d4-tree-view-group")
          ?.querySelector(":scope > .d4-tree-view-group-host");
        const groups = postgresHost
          ? Array.from(postgresHost.children).filter((c) =>
              c.classList.contains("d4-tree-view-group"),
            )
          : [];
        const northwindGroup = groups.find(
          (g) =>
            g
              .querySelector(
                ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
              )
              ?.textContent?.trim() === "Northwind",
        );

        // Browse may already have expanded Northwind; expand if not.
        let nHost = northwindGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        if (!nHost || nHost.children.length === 0) {
          const nRow = northwindGroup?.querySelector(
            ":scope > .d4-tree-view-node",
          ) as HTMLElement | null;
          (
            nRow?.querySelector(
              ".d4-tree-view-tri, .d4-tree-view-group-tri",
            ) as HTMLElement | null
          )?.click();
          await new Promise((r) => setTimeout(r, 4500));
          nHost = northwindGroup?.querySelector(
            ":scope > .d4-tree-view-group-host",
          );
        }

        const nChildren = nHost
          ? Array.from(nHost.children).filter((c) =>
              c.classList.contains("d4-tree-view-group"),
            )
          : [];
        const schemasGroup = nChildren.find(
          (g) =>
            g
              .querySelector(
                ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
              )
              ?.textContent?.trim() === "Schemas",
        );
        let schemasHost = schemasGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        if (!schemasHost || schemasHost.children.length === 0) {
          const sRow = schemasGroup?.querySelector(
            ":scope > .d4-tree-view-node",
          ) as HTMLElement | null;
          sRow?.scrollIntoView({ block: "center" });
          (
            sRow?.querySelector(
              ".d4-tree-view-tri, .d4-tree-view-group-tri",
            ) as HTMLElement | null
          )?.click();
          await new Promise((r) => setTimeout(r, 3500));
          schemasHost = schemasGroup?.querySelector(
            ":scope > .d4-tree-view-group-host",
          );
        }

        const sChildren = schemasHost
          ? Array.from(schemasHost.children).filter((c) =>
              c.classList.contains("d4-tree-view-group"),
            )
          : [];
        const pubGroup = sChildren.find(
          (g) =>
            g
              .querySelector(
                ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
              )
              ?.textContent?.trim() === "public",
        );
        let pubHost = pubGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        if (!pubHost || pubHost.children.length === 0) {
          const pRow = pubGroup?.querySelector(
            ":scope > .d4-tree-view-node",
          ) as HTMLElement | null;
          pRow?.scrollIntoView({ block: "center" });
          (
            pRow?.querySelector(
              ".d4-tree-view-tri, .d4-tree-view-group-tri",
            ) as HTMLElement | null
          )?.click();
          await new Promise((r) => setTimeout(r, 4500));
          pubHost = pubGroup?.querySelector(
            ":scope > .d4-tree-view-group-host",
          );
        }

        const pubChildren = pubHost
          ? Array.from(pubHost.children).filter((c) =>
              c.classList.contains("d4-tree-view-group"),
            )
          : [];
        const customersGroup = pubChildren.find(
          (g) =>
            g
              .querySelector(
                ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
              )
              ?.textContent?.trim() === "customers",
        );
        const cRow = customersGroup?.querySelector(
          ":scope > .d4-tree-view-node",
        ) as HTMLElement | null;
        cRow?.scrollIntoView({ block: "center" });
        cRow?.dispatchEvent(
          new MouseEvent("contextmenu", {
            bubbles: true,
            cancelable: true,
            button: 2,
          }),
        );
        await new Promise((r) => setTimeout(r, 700));
        const items = Array.from(
          document.querySelectorAll(".d4-menu-item-label"),
        )
          .map((el) => el.textContent?.trim())
          .filter((t): t is string => !!t);
        document.body.click();
        await new Promise((r) => setTimeout(r, 200));
        return { customersFound: !!customersGroup, items };
      });
      expect(result.customersFound).toBe(true);
      expect(result.items).toContain("Get All");
      expect(result.items).toContain("Get Top 100");
      expect(result.items).toContain("New SQL Query...");
      expect(result.items).toContain("New Visual Query...");
    },
  );

  if (stepErrors.length > 0) {
    const summary = stepErrors
      .map((e) => `  - ${e.step}: ${e.error}`)
      .join("\n");
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
