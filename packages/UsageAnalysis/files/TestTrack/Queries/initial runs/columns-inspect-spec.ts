import {expect} from '@playwright/test';
import {test} from '../../shared-page';
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

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
        const poll = async (fn: () => boolean, capMs: number) => {
          const deadline = Date.now() + capMs;
          while (Date.now() < deadline) {
            if (fn()) return true;
            await wait(25);
          }
          return false;
        };

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

        const hostOf = (label: Element | null | undefined) =>
          label
            ?.closest(".d4-tree-view-group")
            ?.querySelector(":scope > .d4-tree-view-group-host");

        const stages: string[] = [];
        const tables: {
          table?: string;
          columnCount: number;
          clicks: number;
        }[] = [];

        const waitForSelector = async (sel: string, capMs = 12_000) => {
          await poll(() => !!document.querySelector(sel), capMs);
          return document.querySelector(sel);
        };

        const dbGroup = document.querySelector('[name="tree-Databases"]');
        expandIfNeeded(dbGroup);
        await poll(
          () => !!document.querySelector(`[name="tree-Databases---${provider}"]`),
          700,
        );
        stages.push("Databases");

        const providerNode = await waitForSelector(
          `[name="tree-Databases---${provider}"]`,
        );
        if (!providerNode) return { ok: false, stage: provider, stages };
        expandIfNeeded(providerNode);
        stages.push(provider);

        const connNode = await waitForSelector(
          `[name="tree-Databases---${provider}---${connection}"]`,
        );
        if (!connNode) return { ok: false, stage: connection, stages };
        expandIfNeeded(connNode);
        const connGroup = connNode.closest(".d4-tree-view-group");
        const connHostAt = () =>
          connGroup?.querySelector(":scope > .d4-tree-view-group-host");
        await poll(() => !!findLabelInHost(connHostAt(), "Schemas"), 2500);
        stages.push(connection);

        const connHost = connHostAt();
        const schemasLabel = findLabelInHost(connHost, "Schemas");
        if (!schemasLabel) return { ok: false, stage: "Schemas", stages };
        stages.push("Schemas");
        expandIfNeeded(schemasLabel.closest(".d4-tree-view-node"));
        await poll(
          () => !!findLabelInHost(hostOf(schemasLabel), "public"),
          2500,
        );

        const schemasHost = hostOf(schemasLabel);
        const publicLabel = findLabelInHost(schemasHost, "public");
        if (!publicLabel) return { ok: false, stage: "public", stages };
        stages.push("public");
        expandIfNeeded(publicLabel.closest(".d4-tree-view-node"));
        await poll(() => {
          const h = hostOf(publicLabel);
          return (
            !!h &&
            h.querySelectorAll(":scope > .d4-tree-view-group").length > 0
          );
        }, 3500);

        const publicHost = hostOf(publicLabel);
        const tableGroups = publicHost
          ? Array.from(
              publicHost.querySelectorAll(":scope > .d4-tree-view-group"),
            )
          : [];
        stages.push(`tables:${tableGroups.length}`);

        const w: any = window;
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
          const colsAt = () => {
            const thost = tg.querySelector(
              ":scope > .d4-tree-view-group-host",
            );
            return thost
              ? Array.from(
                  thost.querySelectorAll(
                    ".d4-tree-view-item-label, .d4-tree-view-group-label",
                  ),
                )
              : [];
          };
          await poll(() => colsAt().length > 0, 900);
          const colEls = colsAt();
          let clicks = 0;
          for (const c of colEls) {
            const expected = (c.textContent ?? "").trim();
            (c as HTMLElement).click();
            clicks++;
            // the click sets the current object; wait for that rather than a flat 120ms
            await poll(() => {
              const o = w.grok?.shell?.o;
              return !!o && (o.name ?? String(o)) === expected;
            }, 120);
          }
          tables.push({ table: label, columnCount: colEls.length, clicks });
          if (tri) tri.click();
          await poll(() => colsAt().length === 0, 150);
        }

        // give a balloon raised by the last click its chance to land before reading
        await poll(
          () =>
            Array.from(
              document.querySelectorAll(".d4-balloon, .grok-balloon"),
            ).some((b) => (b as HTMLElement).offsetParent !== null),
          400,
        );
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

  await page.evaluate(async () => {
    const node = document.querySelector(
      '[name="tree-Databases---PostgresDart---NorthwindTest"]',
    ) as HTMLElement | null;
    const tri = node?.querySelector(
      ".d4-tree-view-tri.d4-tree-view-tri-expanded",
    ) as HTMLElement | null;
    if (tri) tri.click();
    const group = node?.closest(".d4-tree-view-group");
    for (let i = 0; i < 20; i++) {
      const host = group?.querySelector(":scope > .d4-tree-view-group-host");
      if (!host || host.children.length === 0) break;
      await new Promise((r) => setTimeout(r, 25));
    }
  });

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
