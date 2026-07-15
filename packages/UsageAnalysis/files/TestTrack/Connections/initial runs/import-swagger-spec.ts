import { test, expect } from "@playwright/test";
import * as fs from "fs";
import * as path from "path";
import {
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

const SWAGGER_PATH = path.resolve(
  __dirname,
  "../../../../Samples/swaggers/openweathermap.yaml",
);

test("Connections / Import SWAGGER", async ({ page }) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const existing = await (window as any).grok.dapi.connections
      .filter('name = "OpenWeatherMap"')
      .list();
    for (const c of existing) {
      try {
        await (window as any).grok.dapi.connections.delete(c);
      } catch {}
    }
  });

  await softStep(
    "Step 1-2: locate and read openweathermap.yaml from Samples package",
    async () => {
      const exists = fs.existsSync(SWAGGER_PATH);
      expect(exists, `swagger file expected at ${SWAGGER_PATH}`).toBe(true);
      const content = fs.readFileSync(SWAGGER_PATH, "utf-8");
      expect(content).toContain("OpenWeatherMap");
      (test.info() as any).swaggerContent = content;
    },
  );

  await softStep("Step 3: drag-and-drop yaml onto Datagrok", async () => {
    const yamlText = fs.readFileSync(SWAGGER_PATH, "utf-8");
    const result = await page.evaluate(async (content: string) => {
      const fsApi =
        (window as any).webkitRequestFileSystem ||
        (window as any).requestFileSystem;
      if (!fsApi)
        return { ok: false, reason: "webkitRequestFileSystem unavailable" };
      const realEntry: any = await new Promise((resolve, reject) => {
        fsApi.call(
          window,
          (window as any).TEMPORARY,
          1024 * 1024,
          (fsys: any) => {
            fsys.root.getFile(
              "openweathermap.yaml",
              { create: true },
              (entry: any) => {
                entry.createWriter((writer: any) => {
                  writer.onwriteend = () => {
                    entry.createWriter((w: any) => {
                      w.onwriteend = () => resolve(entry);
                      w.onerror = (e: any) => reject(e);
                      w.write(
                        new Blob([content], { type: "application/x-yaml" }),
                      );
                    }, reject);
                  };
                  writer.onerror = (e: any) => reject(e);
                  writer.truncate(0);
                }, reject);
              },
              reject,
            );
          },
          reject,
        );
      });
      const file: File = await new Promise((resolve, reject) =>
        realEntry.file(resolve, reject),
      );
      const dt = new DataTransfer();
      dt.items.add(file);
      const orig = (DataTransferItem.prototype as any).webkitGetAsEntry;
      (DataTransferItem.prototype as any).webkitGetAsEntry = function () {
        return this.kind === "file" ? realEntry : null;
      };
      try {
        const xpRoot = document.getElementById("rootDiv");
        if (!xpRoot) return { ok: false, reason: "rootDiv not found" };
        xpRoot.dispatchEvent(
          new DragEvent("dragenter", {
            bubbles: true,
            cancelable: true,
            dataTransfer: dt,
          }),
        );
        await new Promise((r) => setTimeout(r, 400));
        const overlay: any = Array.from(xpRoot.children).find(
          (c: any) =>
            c.tagName === "DIV" && c.style && c.style.zIndex === "10000",
        );
        if (!overlay) return { ok: false, reason: "overlay not created" };
        overlay.dispatchEvent(
          new DragEvent("dragover", {
            bubbles: true,
            cancelable: true,
            dataTransfer: dt,
          }),
        );
        overlay.dispatchEvent(
          new DragEvent("drop", {
            bubbles: true,
            cancelable: true,
            dataTransfer: dt,
          }),
        );
        await new Promise((r) => setTimeout(r, 8000));
      } finally {
        (DataTransferItem.prototype as any).webkitGetAsEntry = orig;
      }
      const conns = await (window as any).grok.dapi.connections
        .filter('name = "OpenWeatherMap"')
        .list();
      return { ok: conns.length > 0, count: conns.length };
    }, yamlText);
    expect(
      result.ok,
      `swagger import expected to create OpenWeatherMap connection (${(result as any).reason ?? ""})`,
    ).toBe(true);
  });

  await softStep(
    "Step 4: navigate Browse > Platform > Functions > OpenAPI > OpenWeatherMap",
    async () => {
      await page.goto(process.env.DATAGROK_URL + "/browse", {
        waitUntil: "networkidle",
        timeout: 30_000,
      });
      const found = await page.evaluate(async () => {
        const expand = async (labelText: string) => {
          for (let i = 0; i < 60; i++) {
            const label = Array.from(
              document.querySelectorAll(".d4-tree-view-group-label"),
            ).find((l) => l.textContent?.trim() === labelText) as
              | HTMLElement
              | undefined;
            if (label) {
              const row = label.closest(".d4-tree-view-node");
              row?.scrollIntoView({ block: "center" });
              const tri = row?.querySelector(
                ".d4-tree-view-tri, .d4-tree-view-group-tri",
              ) as HTMLElement | null;
              tri?.click();
              return true;
            }
            await new Promise((r) => setTimeout(r, 250));
          }
          return false;
        };
        if (!(await expand("Platform"))) return { stage: "no Platform" };
        await new Promise((r) => setTimeout(r, 1500));
        if (!(await expand("Functions"))) return { stage: "no Functions" };
        await new Promise((r) => setTimeout(r, 1500));
        if (!(await expand("OpenAPI"))) return { stage: "no OpenAPI" };
        await new Promise((r) => setTimeout(r, 3000));
        for (let i = 0; i < 25; i++) {
          const has = Array.from(
            document.querySelectorAll(".d4-tree-view-group-label"),
          ).some((l) => l.textContent?.trim() === "OpenWeatherMap");
          if (has) return { stage: "found" };
          await new Promise((r) => setTimeout(r, 300));
        }
        return { stage: "OpenWeatherMap not visible" };
      });
      expect(
        found.stage,
        `Browse > Platform > Functions > OpenAPI: ${JSON.stringify(found)}`,
      ).toBe("found");
    },
  );

  await softStep("Step 5: right-click connection, select Edit", async () => {
    const result = await page.evaluate(async () => {
      const target = Array.from(
        document.querySelectorAll(".d4-tree-view-group-label"),
      ).find((l) => l.textContent?.trim() === "OpenWeatherMap") as
        | HTMLElement
        | undefined;
      if (!target) return { stage: "no target" };
      const row = target.closest(".d4-tree-view-node") as HTMLElement | null;
      row?.scrollIntoView({ block: "center" });
      const r = row?.getBoundingClientRect();
      row?.dispatchEvent(
        new MouseEvent("contextmenu", {
          bubbles: true,
          cancelable: true,
          button: 2,
          clientX: (r?.left ?? 0) + 10,
          clientY: (r?.top ?? 0) + 5,
        }),
      );
      await new Promise((rr) => setTimeout(rr, 700));
      const editItem = (document.querySelector('[name="div-Edit..."]') ??
        Array.from(document.querySelectorAll(".d4-menu-item-label")).find(
          (el) => el.textContent?.trim() === "Edit...",
        )) as HTMLElement | null;
      if (!editItem) return { stage: "no Edit menu item" };
      editItem.click();
      for (let i = 0; i < 50; i++) {
        await new Promise((rr) => setTimeout(rr, 200));
        const dlg = document.querySelector(".d4-dialog");
        if (dlg && dlg.querySelector('[name="input-ApiKey"]'))
          return { stage: "ready" };
      }
      const dlg = document.querySelector(".d4-dialog");
      const title = dlg
        ?.querySelector(".d4-dialog-title, .d4-dialog-header")
        ?.textContent?.trim();
      const inputs = Array.from(dlg?.querySelectorAll("input[name]") ?? []).map(
        (i: any) => i.name,
      );
      return {
        stage: "no apikey input ready",
        dialogPresent: !!dlg,
        title,
        inputs,
      };
    });
    expect(
      result.stage,
      `Edit Connection dialog: ${JSON.stringify(result)}`,
    ).toBe("ready");
  });

  await softStep(
    "Step 6: enter ApiKey (placeholder for QA-supplied key)",
    async () => {
      const result = await page.evaluate(async () => {
        const dialog = document.querySelector(".d4-dialog");
        if (!dialog) return { stage: "no dialog" };
        const apiKeyInput: HTMLInputElement | null = dialog.querySelector(
          '[name="input-ApiKey"]',
        );
        if (!apiKeyInput)
          return {
            stage: "no apikey input",
            labels: Array.from(dialog.querySelectorAll("label"))
              .map((l) => l.textContent?.trim())
              .filter(Boolean),
          };
        const ns = Object.getOwnPropertyDescriptor(
          HTMLInputElement.prototype,
          "value",
        )!.set!;
        ns.call(apiKeyInput, "PLACEHOLDER_NOT_A_REAL_KEY");
        apiKeyInput.dispatchEvent(new Event("input", { bubbles: true }));
        apiKeyInput.dispatchEvent(new Event("change", { bubbles: true }));
        const okBtn: any =
          dialog.querySelector('[name="button-OK"]') ??
          document.querySelector('[name="button-OK"]');
        if (!okBtn) return { stage: "no OK button" };
        okBtn.click();
        for (let i = 0; i < 60; i++) {
          await new Promise((r) => setTimeout(r, 300));
          if (!document.querySelector(".d4-dialog"))
            return { stage: "closed", i };
        }
        return { stage: "still open after 18s" };
      });
      expect(result.stage, `Edit dialog flow: ${JSON.stringify(result)}`).toBe(
        "closed",
      );
    },
  );

  await softStep("Step 7: run all queries on OpenWeatherMap", async () => {
    const results = await page.evaluate(async () => {
      const conn = (
        await (window as any).grok.dapi.connections
          .filter('name = "OpenWeatherMap"')
          .list()
      )[0];
      const queries = await (window as any).grok.dapi.queries
        .filter('connection.id = "' + conn.id + '"')
        .list();
      const out: any[] = [];
      for (const q of queries) {
        const params: any = {};
        for (const p of q.inputs || []) {
          if (p.name === "q") params.q = "London";
          else if (p.name === "lat") params.lat = 51.5;
          else if (p.name === "lon") params.lon = -0.1;
          else if (p.name === "cnt") params.cnt = 1;
          else if (p.name === "bbox") params.bbox = [-1, 50, 1, 52, 5];
          else if (p.name === "start") params.start = "2024-01-01T00:00:00Z";
          else if (p.name === "end") params.end = "2024-01-02T00:00:00Z";
        }
        try {
          await q
            .prepare(params)
            .call(true, null, { processed: false, report: false });
          out.push({ name: q.friendlyName, ok: true });
        } catch (e: any) {
          out.push({
            name: q.friendlyName,
            ok: false,
            error: String(e?.message ?? e).split("\n")[0],
          });
        }
      }
      return out;
    });
    expect(
      results.length,
      "expected to invoke all 7 OpenAPI-derived queries",
    ).toBeGreaterThanOrEqual(7);
    for (const r of results)
      expect(
        r.ok,
        `query ${r.name} expected to be invocable: ${(r as any).error ?? ""}`,
      ).toBe(true);
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors
      .map((e) => `  - ${e.step}: ${e.error}`)
      .join("\n");
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
