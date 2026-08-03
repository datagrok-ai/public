import { test, expect } from "@playwright/test";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

test("Queries: browse and save project", async ({ page }) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const w = window as any;
    document.body.classList.add("selenium");
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
  });

  await softStep("Step 1: Navigate to Databases > Postgres", async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      // Verify both CHEMBL and Northwind connections exist on Postgres provider
      const conns = await w.grok.dapi.connections.list({ pageSize: 1000 });
      const pg = conns.filter((c: any) => c.dataSource === "Postgres");
      const chembl = pg.find(
        (c: any) => c.friendlyName === "CHEMBL" || c.name === "Chembl",
      );
      const nw = pg.find(
        (c: any) =>
          c.friendlyName === "Northwind" || c.name === "PostgresNorthwind",
      );
      return { chembl: !!chembl, nw: !!nw };
    });
    expect(result.chembl).toBe(true);
    expect(result.nw).toBe(true);
  });

  await softStep(
    "Step 2: Preview + run sample queries on Northwind and CHEMBL",
    async () => {
      const result = await page.evaluate(async () => {
        const w = window as any;
        const conns = await w.grok.dapi.connections.list({ pageSize: 1000 });
        const nwConn = conns.find((c: any) => c.name === "PostgresNorthwind");
        const chConn = conns.find(
          (c: any) => c.name === "Chembl" && c.dataSource === "Postgres",
        );
        const all = await w.grok.dapi.queries.list({ pageSize: 1000 });
        const nwQ = all.find(
          (x: any) =>
            x.connection?.id === nwConn?.id && x.name === "PostgresProducts",
        );
        const chQ = all.find(
          (x: any) =>
            x.connection?.id === chConn?.id && (x.inputs || []).length === 0,
        );
        const nwDf = await nwQ.executeTable();
        const chDf = await chQ.executeTable();
        w.grok.shell.addTableView(nwDf);
        await new Promise((r) => setTimeout(r, 600));
        w.grok.shell.addTableView(chDf);
        return { nwRows: nwDf.rowCount, chRows: chDf.rowCount };
      });
      expect(result.nwRows).toBeGreaterThan(0);
      expect(result.chRows).toBeGreaterThan(0);
    },
  );

  await softStep(
    "Step 3: Open FRAC query, verify Context Panel tabs",
    async () => {
      const tabs = await page.evaluate(async () => {
        const w = window as any;
        w.grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 800));
        const all = await w.grok.dapi.queries.list({ pageSize: 1000 });
        const q = all.find(
          (x: any) => x.name === "FracClassificationWithSubstructure",
        );
        // Select the query (browse-style preview)
        w.grok.shell.o = q;
        await new Promise((r) => setTimeout(r, 1500));
        const headers = Array.from(
          document.querySelectorAll(".d4-accordion-pane-header"),
        ).map((h) => h.textContent?.trim() || "");
        return headers;
      });
      expect(tabs.some((t) => t.startsWith("Run"))).toBe(true);
      expect(tabs.some((t) => t.startsWith("Query"))).toBe(true);
      expect(tabs.some((t) => t.startsWith("Transformations"))).toBe(true);
    },
  );

  await softStep(
    "Step 4: Open FRAC run dialog; first param change clears dependent params",
    async () => {
      const result = await page.evaluate(async () => {
        const runHdr = Array.from(
          document.querySelectorAll(".d4-accordion-pane-header"),
        ).find((h) => h.textContent?.trim().toLowerCase().startsWith("run")) as
          | HTMLElement
          | undefined;
        runHdr?.click();
        await new Promise((r) => setTimeout(r, 800));
        const runBtn = document.querySelector<HTMLElement>(
          '[name="button-RUN"]',
        );
        runBtn?.click();
        await new Promise((r) => setTimeout(r, 2000));
        const dlg = document.querySelector(".d4-dialog");
        const sel1 = dlg?.querySelector<HTMLSelectElement>(
          '[name="input-host-Level1"] select',
        );
        const opts = Array.from(sel1?.options || []).map((o) => o.value);
        const newVal = opts.find((v) => v && v !== sel1?.value);
        if (sel1 && newVal) {
          sel1.value = newVal;
          sel1.dispatchEvent(new Event("change", { bubbles: true }));
        }
        await new Promise((r) => setTimeout(r, 1500));
        const vals = ["Level1", "Level2", "Level3", "Level4"].map((n) => {
          const h = dlg?.querySelector(
            `[name="input-host-${n}"] select`,
          ) as HTMLSelectElement | null;
          return { n, v: h?.value };
        });
        return { hasDialog: !!dlg, vals };
      });
      expect(result.hasDialog).toBe(true);
      expect(result.vals.find((v) => v.n === "Level2")?.v).toBe("");
      expect(result.vals.find((v) => v.n === "Level3")?.v).toBe("");
      expect(result.vals.find((v) => v.n === "Level4")?.v).toBe("");
    },
  );

  await softStep(
    "Step 5: Set parameters and click OK; result table opens",
    async () => {
      const result = await page.evaluate(async () => {
        const w = window as any;
        const dlg = document.querySelector(".d4-dialog");
        const sel1 = dlg?.querySelector<HTMLSelectElement>(
          '[name="input-host-Level1"] select',
        );
        if (sel1) {
          sel1.value = "STEROL BIOSYNTHESIS IN MEMBRANES";
          sel1.dispatchEvent(new Event("change", { bubbles: true }));
        }
        await new Promise((r) => setTimeout(r, 1500));
        const sel2 = dlg?.querySelector<HTMLSelectElement>(
          '[name="input-host-Level2"] select',
        );
        if (sel2) {
          const opts = Array.from(sel2.options)
            .map((o) => o.value)
            .filter(Boolean);
          if (opts.length) {
            sel2.value = opts[0];
            sel2.dispatchEvent(new Event("change", { bubbles: true }));
          }
        }
        await new Promise((r) => setTimeout(r, 1500));
        const ok = dlg?.querySelector<HTMLElement>('[name="button-OK"]');
        ok?.click();
        for (let i = 0; i < 30; i++) {
          await new Promise((r) => setTimeout(r, 1000));
          if (!document.querySelector(".d4-dialog")) break;
        }
        await new Promise((r) => setTimeout(r, 3000));
        const tv = w.grok.shell.tv;
        return {
          viewType: w.grok.shell.v?.type,
          rows: tv?.dataFrame?.rowCount,
        };
      });
      expect(result.viewType).toBe("TableView");
      expect(result.rows).toBeGreaterThan(0);
    },
  );

  await softStep("Step 6: Add Trellis plot to result view", async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const tp = document.querySelector<HTMLElement>(
        '[name="icon-trellis-plot"]',
      );
      tp?.click();
      await new Promise((r) => setTimeout(r, 2500));
      const tv = w.grok.shell.tv;
      const viewers =
        tv && tv.viewers ? Array.from(tv.viewers).map((v: any) => v.type) : [];
      return { viewers };
    });
    expect(result.viewers).toContain("Trellis plot");
  });

  await softStep("Step 7: Save project (with embedded layout)", async () => {
    const result = await page.evaluate(async () => {
      const w = window as any;
      const saveBtn = document.querySelector<HTMLElement>(
        '[name="button-Save"]',
      );
      saveBtn?.click();
      await new Promise((r) => setTimeout(r, 3000));
      const dlg = document.querySelector(".d4-dialog");
      const inputs = Array.from(
        dlg?.querySelectorAll('input[type="text"], input:not([type])') || [],
      );
      const nameInput = inputs.find((i) =>
        /FRAC classification with substructure/i.test(
          (i as HTMLInputElement).value || "",
        ),
      ) as HTMLInputElement | undefined;
      const newName = "FRAC scenario test " + Date.now();
      if (nameInput) {
        nameInput.focus();
        const proto = Object.getPrototypeOf(nameInput);
        const setter = Object.getOwnPropertyDescriptor(proto, "value")?.set;
        setter?.call(nameInput, newName);
        nameInput.dispatchEvent(new Event("input", { bubbles: true }));
        nameInput.dispatchEvent(new Event("change", { bubbles: true }));
      }
      await new Promise((r) => setTimeout(r, 600));
      const ok = dlg?.querySelector<HTMLElement>('[name="button-OK"]');
      ok?.click();
      for (let i = 0; i < 30; i++) {
        await new Promise((r) => setTimeout(r, 1000));
        if (!document.querySelector(".d4-dialog")) break;
      }
      await new Promise((r) => setTimeout(r, 3000));
      const proj = await w.grok.dapi.projects
        .filter("FRAC scenario test")
        .first();
      return { savedId: proj?.id, savedName: proj?.name };
    });
    expect(result.savedId).toBeTruthy();
  });

  await softStep(
    "Step 8: Close all and reopen project; trellis plot persists",
    async () => {
      const result = await page.evaluate(async () => {
        const w = window as any;
        w.grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 2000));
        const proj = await w.grok.dapi.projects
          .filter("FRAC scenario test")
          .first();
        if (!proj) return { viewers: [], rows: 0 };
        await proj.open();
        for (let i = 0; i < 60; i++) {
          await new Promise((r) => setTimeout(r, 1000));
          const tv = w.grok.shell.tv;
          if (tv && tv.viewers && Array.from(tv.viewers).length > 1) break;
        }
        await new Promise((r) => setTimeout(r, 4000));
        const tv = w.grok.shell.tv;
        const viewers =
          tv && tv.viewers
            ? Array.from(tv.viewers).map((v: any) => v.type)
            : [];
        return { viewers, rows: tv?.dataFrame?.rowCount };
      });
      expect(result.viewers).toContain("Trellis plot");
      expect(result.viewers).toContain("Grid");
      expect(result.rows).toBeGreaterThan(0);
    },
  );

  // Cleanup
  await page.evaluate(async () => {
    const w = window as any;
    const old = await w.grok.dapi.projects.filter("FRAC scenario test").first();
    if (old) await w.grok.dapi.projects.delete(old);
    w.grok.shell.closeAll();
  });

  if (stepErrors.length > 0)
    throw new Error(
      "Step failures:\n" +
        stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join("\n"),
    );
});
