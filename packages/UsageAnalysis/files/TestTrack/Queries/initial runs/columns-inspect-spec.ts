import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

/**
 * Walks Browse → Databases → {Provider} → {Connection} → Schemas → public,
 * expands every public table, clicks every column, and asserts that no error
 * balloons appear and the Context Panel accordion has no error indicators.
 *
 * The Schemas/public/table group nodes have no `name=` attributes, so they
 * must be located by label text scoped under the parent's children host.
 */
test("Queries — columns inspection on PostgresDart NorthwindTest and Postgres Northwind", async ({
  page,
}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });
  await page.locator('[name="Browse"]').waitFor({ timeout: 30_000 });

  await page.goto(`${process.env.DATAGROK_URL}/browse`);
  await page.locator(".d4-tree-view-root").waitFor({ timeout: 30_000 });

  const walkAndInspect = async (provider: string, connection: string) => {
    return await page.evaluate(
      async ({
        provider,
        connection,
      }: {
        provider: string;
        connection: string;
      }) => {
        const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));

        const expandIfNeeded = (node: Element | null) => {
          if (!node) return false;
          const tri = node.querySelector(
            ":scope > .d4-tree-view-tri",
          ) as HTMLElement | null;
          if (tri && !tri.classList.contains("d4-tree-view-tri-expanded")) {
            tri.click();
            return true;
          }
          if (!tri)
            (
              node.querySelector(
                ".d4-tree-view-group-label, .d4-tree-view-item-label",
              ) as HTMLElement
            )?.click();
          return false;
        };

        const findLabelInHost = (
          host: Element | null | undefined,
          text: string,
        ) =>
          host &&
          (Array.from(
            host.querySelectorAll(
              ":scope > .d4-tree-view-group > .d4-tree-view-node > .d4-tree-view-group-label",
            ),
          ).find((e) => e.textContent?.trim() === text) as
            | HTMLElement
            | undefined);

        const stages: string[] = [];
        const tables: {
          table?: string;
          columnCount: number;
          clicks: number;
        }[] = [];

        const waitForSelector = async (sel: string, maxIter = 40) => {
          for (let i = 0; i < maxIter; i++) {
            const el = document.querySelector(sel);
            if (el) return el;
            await wait(300);
          }
          return null;
        };

        // Make sure Databases group is expanded.
        const dbGroup = document.querySelector('[name="tree-Databases"]');
        expandIfNeeded(dbGroup);
        await wait(700);
        stages.push("Databases");

        // Expand provider node and wait for the connection child to materialize.
        const providerNode = await waitForSelector(
          `[name="tree-Databases---${provider}"]`,
        );
        if (!providerNode) return { ok: false, stage: provider, stages };
        expandIfNeeded(providerNode);
        stages.push(provider);

        // The connection child appears asynchronously after provider expands.
        const connNode = await waitForSelector(
          `[name="tree-Databases---${provider}---${connection}"]`,
        );
        if (!connNode) return { ok: false, stage: connection, stages };
        expandIfNeeded(connNode);
        await wait(2500);
        stages.push(connection);

        // Find Schemas under the connection's host.
        const connGroup = connNode.closest(".d4-tree-view-group");
        const connHost = connGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        const schemasLabel = findLabelInHost(connHost, "Schemas");
        if (!schemasLabel) return { ok: false, stage: "Schemas", stages };
        stages.push("Schemas");
        expandIfNeeded(schemasLabel.closest(".d4-tree-view-node"));
        await wait(2500);

        // Find public under Schemas' host.
        const schemasGroup = schemasLabel.closest(".d4-tree-view-group");
        const schemasHost = schemasGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        const publicLabel = findLabelInHost(schemasHost, "public");
        if (!publicLabel) return { ok: false, stage: "public", stages };
        stages.push("public");
        expandIfNeeded(publicLabel.closest(".d4-tree-view-node"));
        await wait(3500);

        const publicGroup = publicLabel.closest(".d4-tree-view-group");
        const publicHost = publicGroup?.querySelector(
          ":scope > .d4-tree-view-group-host",
        );
        const tableGroups = publicHost
          ? Array.from(
              publicHost.querySelectorAll(":scope > .d4-tree-view-group"),
            )
          : [];
        stages.push(`tables:${tableGroups.length}`);

        // Expand each table, click each column, then collapse.
        for (const tg of tableGroups) {
          const label = tg
            .querySelector(
              ":scope > .d4-tree-view-node > .d4-tree-view-group-label",
            )
            ?.textContent?.trim();
          const tnode = tg.querySelector(":scope > .d4-tree-view-node");
          const tri = tnode?.querySelector(
            ".d4-tree-view-tri",
          ) as HTMLElement | null;
          if (tri) tri.click();
          else (tnode as HTMLElement | null)?.click();
          await wait(900);
          const thost = tg.querySelector(":scope > .d4-tree-view-group-host");
          const colEls = thost
            ? Array.from(
                thost.querySelectorAll(
                  ".d4-tree-view-item-label, .d4-tree-view-group-label",
                ),
              )
            : [];
          let clicks = 0;
          for (const c of colEls) {
            (c as HTMLElement).click();
            clicks++;
            await wait(120);
          }
          tables.push({ table: label, columnCount: colEls.length, clicks });
          if (tri) tri.click();
          await wait(150);
        }

        const balloons = Array.from(
          document.querySelectorAll(".d4-balloon, .grok-balloon"),
        )
          .filter((b) => (b as HTMLElement).offsetParent !== null)
          .map((b) => (b.textContent?.trim() ?? "").slice(0, 160));
        const accordionErrors = Array.from(
          document.querySelectorAll(".grok-prop-panel .d4-accordion-pane"),
        )
          .filter(
            (a) =>
              !!a.querySelector(".d4-error, .grok-error, .d4-accordion-error"),
          )
          .map((a) =>
            a.querySelector(".d4-accordion-pane-header")?.textContent?.trim(),
          );

        return {
          ok: true,
          stages,
          tableCount: tables.length,
          totalClicks: tables.reduce((s, t) => s + t.clicks, 0),
          balloons,
          accordionErrors,
        };
      },
      { provider, connection },
    );
  };

  // ---------- Part 1: PostgresDart → NorthwindTest ----------
  let part1: any;
  await softStep(
    "Browse → Databases → PostgresDart → NorthwindTest → Schemas → public",
    async () => {
      part1 = await walkAndInspect("PostgresDart", "NorthwindTest");
      expect(part1.ok, JSON.stringify(part1)).toBe(true);
      expect(part1.stages).toContain("public");
      expect(part1.tableCount).toBeGreaterThan(0);
    },
  );

  await softStep(
    "Part 1 — expand each DB table, click each column, no errors on Context Panel",
    async () => {
      expect(part1?.totalClicks ?? 0).toBeGreaterThan(0);
      expect(part1?.balloons ?? []).toEqual([]);
      expect(part1?.accordionErrors ?? []).toEqual([]);
    },
  );

  // ---------- Part 2: Postgres → Northwind ----------
  // Collapse Part 1 connection to keep the DOM small.
  await page.evaluate(() => {
    const node = document.querySelector(
      '[name="tree-Databases---PostgresDart---NorthwindTest"]',
    ) as HTMLElement | null;
    const tri = node?.querySelector(
      ".d4-tree-view-tri.d4-tree-view-tri-expanded",
    ) as HTMLElement | null;
    if (tri) tri.click();
  });
  await page.waitForTimeout(500);

  let part2: any;
  await softStep(
    "Browse → Databases → Postgres → Northwind → Schemas → public",
    async () => {
      part2 = await walkAndInspect("Postgres", "Northwind");
      expect(part2.ok, JSON.stringify(part2)).toBe(true);
      expect(part2.stages).toContain("public");
      expect(part2.tableCount).toBeGreaterThan(0);
    },
  );

  await softStep(
    "Part 2 — expand each DB table, click each column, no errors on Context Panel",
    async () => {
      expect(part2?.totalClicks ?? 0).toBeGreaterThan(0);
      expect(part2?.balloons ?? []).toEqual([]);
      expect(part2?.accordionErrors ?? []).toEqual([]);
    },
  );

  if (stepErrors.length > 0)
    throw new Error(
      "Soft step failures:\n" +
        stepErrors.map((e) => `- ${e.step}: ${e.error}`).join("\n"),
    );
});
