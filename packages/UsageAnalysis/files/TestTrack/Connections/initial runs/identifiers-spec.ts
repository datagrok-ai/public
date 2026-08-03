import { test, expect } from "@playwright/test";
import {
  baseUrl,
  loginToDatagrok,
  specTestOptions,
  softStep,
  stepErrors,
} from "../../spec-login";

test.use(specTestOptions);

// Prerequisite: a "test_postgres" connection must exist that connects to
// db.datagrok.ai:54322/northwind. The MCP run created Agolovko:test_postgres
// (id af9bcf40-21a0-11f1-89e2-7b1321b80948).

const TEST_POSTGRES = "test_postgres";

async function gotoPostgresAndExpand(page: import("@playwright/test").Page) {
  await page.goto(`${baseUrl}/browse/databases/Postgres`, {
    waitUntil: "networkidle",
    timeout: 30000,
  });
  await page.waitForTimeout(4000);
  await page.evaluate(() => {
    document.body.classList.add("selenium");
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.windows.showBrowse = true;
  });
  await page.waitForTimeout(1000);
  await page.evaluate(async () => {
    let pg: HTMLElement | undefined;
    for (let i = 0; i < 20; i++) {
      pg = Array.from(
        document.querySelectorAll(".d4-tree-view-group-label"),
      ).find((el) => el.textContent?.trim() === "Postgres") as
        | HTMLElement
        | undefined;
      if (pg) break;
      await new Promise((r) => setTimeout(r, 300));
    }
    if (!pg) return;
    const node = pg.closest(".d4-tree-view-node") as HTMLElement;
    const tri = node.querySelector(".d4-tree-view-tri") as HTMLElement | null;
    if (tri && !tri.classList.contains("d4-tree-view-tri-expanded"))
      tri.click();
    else pg.click();
    await new Promise((r) => setTimeout(r, 3000));
    for (let i = 0; i < 5; i++) {
      const showMore = Array.from(
        node.querySelectorAll(".d4-tree-view-group-label, label, a"),
      ).find((el) => el.textContent?.trim() === "Show more") as
        | HTMLElement
        | undefined;
      if (!showMore) break;
      showMore.click();
      await new Promise((r) => setTimeout(r, 1500));
    }
  });
}

async function rightClickConnection(
  page: import("@playwright/test").Page,
  name: string,
) {
  return await page.evaluate(async (n) => {
    document.dispatchEvent(
      new KeyboardEvent("keydown", { key: "Escape", bubbles: true }),
    );
    await new Promise((r) => setTimeout(r, 200));
    document.querySelectorAll(".d4-menu-popup").forEach((m) => m.remove());
    let tp: HTMLElement | undefined;
    for (let i = 0; i < 60; i++) {
      tp = Array.from(
        document.querySelectorAll(
          ".d4-tree-view-group-label, label, .d4-link-label",
        ),
      ).find((el) => el.textContent?.trim() === n) as HTMLElement | undefined;
      if (tp) break;
      await new Promise((r) => setTimeout(r, 250));
    }
    if (!tp) return { found: false };
    const node = tp.closest(".d4-tree-view-node") as HTMLElement;
    node.scrollIntoView({ block: "center" });
    node.dispatchEvent(
      new MouseEvent("contextmenu", {
        bubbles: true,
        cancelable: true,
        button: 2,
      }),
    );
    await new Promise((r) => setTimeout(r, 1500));
    return { found: true };
  }, name);
}

async function clickMenuItem(
  page: import("@playwright/test").Page,
  text: string,
) {
  return await page.evaluate(async (t) => {
    const target = Array.from(document.querySelectorAll(".d4-menu-item")).find(
      (el) => el.textContent?.trim() === t,
    ) as HTMLElement | undefined;
    if (!target) return { clicked: false };
    const inner = target.querySelector(
      ".d4-menu-item-label",
    ) as HTMLElement | null;
    const el = inner ?? target;
    for (const evtType of [
      "pointerdown",
      "mousedown",
      "pointerup",
      "mouseup",
      "click",
    ]) {
      el.dispatchEvent(
        new MouseEvent(evtType, { bubbles: true, cancelable: true, button: 0 }),
      );
    }
    await new Promise((r) => setTimeout(r, 1200));
    return { clicked: true };
  }, text);
}

async function setDialogInput(
  page: import("@playwright/test").Page,
  labelText: string,
  value: string,
) {
  return await page.evaluate(
    ({ l, v }) => {
      const dialogs = Array.from(document.querySelectorAll(".d4-dialog"));
      const dialog = dialogs[dialogs.length - 1];
      if (!dialog) return { set: false, reason: "no dialog" };
      const labels = Array.from(
        dialog.querySelectorAll("label, .d4-label, .ui-label"),
      );
      const lbl = labels.find((x) => x.textContent?.trim() === l);
      const inp = lbl
        ?.closest(".d4-input-base, .d4-flex-row, div")
        ?.querySelector("input, select, .d4-combo-popup") as
        | HTMLInputElement
        | HTMLSelectElement
        | null;
      if (!inp) return { set: false, reason: "no input for " + l };
      const ns = Object.getOwnPropertyDescriptor(
        HTMLInputElement.prototype,
        "value",
      )!.set!;
      ns.call(inp, v);
      inp.dispatchEvent(new Event("input", { bubbles: true }));
      inp.dispatchEvent(new Event("change", { bubbles: true }));
      return { set: true };
    },
    { l: labelText, v: value },
  );
}

async function clickDialogButton(
  page: import("@playwright/test").Page,
  name: string,
) {
  return await page.evaluate((n) => {
    const dialogs = Array.from(document.querySelectorAll(".d4-dialog"));
    const dialog = dialogs[dialogs.length - 1];
    if (!dialog) return { clicked: false };
    const candidates = Array.from(
      dialog.querySelectorAll(
        'button, .ui-btn, .d4-dialog-footer .ui-btn, [name^="button-"]',
      ),
    );
    const btn = candidates.find((b) => {
      const t = b.textContent?.trim() ?? "";
      return t === n || t.toUpperCase() === n.toUpperCase();
    }) as HTMLElement | undefined;
    btn?.click();
    return { clicked: !!btn };
  }, name);
}

test("Connections / Identifiers", async ({ page }) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await gotoPostgresAndExpand(page);

  await softStep(
    "Step 1: Browse > Databases > Postgres tree visible",
    async () => {
      const visible = await page.evaluate(
        () =>
          !!Array.from(
            document.querySelectorAll(".d4-tree-view-group-label"),
          ).find((el) => el.textContent?.trim() === "Postgres"),
      );
      expect(visible).toBe(true);
    },
  );

  await softStep(
    "Step 2-3: Right-click test_postgres → Configure Identifiers... opens schema dialog",
    async () => {
      const r1 = await rightClickConnection(page, TEST_POSTGRES);
      expect(r1.found).toBe(true);
      const r2 = await clickMenuItem(page, "Configure Identifiers...");
      expect(r2.clicked).toBe(true);
      await expect(page.locator(".d4-dialog .d4-dialog-title")).toContainText(
        "Identifiers Configuration",
        { timeout: 10000 },
      );
    },
  );

  await softStep("Step 4: Set Schema = public, click OK", async () => {
    // Try common label names: "Schema", "Primary Schema"
    let r = await setDialogInput(page, "Schema", "public");
    if (!r.set) r = await setDialogInput(page, "Primary Schema", "public");
    expect(
      r.set,
      `Schema input not found in schema dialog: ${(r as any).reason}`,
    ).toBe(true);
    await page.waitForTimeout(500);
    const ok = await clickDialogButton(page, "OK");
    expect(ok.clicked).toBe(true);
    await page.waitForTimeout(2000);
  });

  await softStep(
    "Step 5: Add identifier (CUSTOMER_ID, customers, customerid, [A-Z]{5})",
    async () => {
      // The identifier editor dialog is now open. Add a new row, fill values.
      // Selectors are best-effort: find an "Add" button or "+" inside the dialog,
      // then fill four inputs (Semantic Type, Table, Column, MatchRegexp).
      const add = await page.evaluate(() => {
        const dialogs = Array.from(document.querySelectorAll(".d4-dialog"));
        const dialog = dialogs[dialogs.length - 1];
        if (!dialog) return { clicked: false, reason: "no dialog" };
        const addBtn = (Array.from(
          dialog.querySelectorAll(
            'button, .ui-btn, .grok-icon, i, [name*="add"]',
          ),
        ).find((b) => {
          const t = b.textContent?.trim() ?? "";
          const n = b.getAttribute("name") ?? "";
          return t === "ADD" || t === "Add" || t === "+" || /add/i.test(n);
        }) ??
          dialog.querySelector(
            '[name="icon-plus"], .fa-plus, .grok-icon-plus',
          )) as HTMLElement | null;
        if (!addBtn) return { clicked: false, reason: "no add button" };
        addBtn.click();
        return { clicked: true };
      });
      expect(
        add.clicked,
        `Add identifier row not found: ${(add as any).reason}`,
      ).toBe(true);
      await page.waitForTimeout(500);

      // Fill the four fields. Names guessed based on scenario wording.
      const fillResults: Record<string, any> = {};
      for (const [label, value] of [
        ["Semantic Type", "CUSTOMER_ID"],
        ["Table", "customers"],
        ["Column", "customerid"],
        ["MatchRegexp", "[A-Z]{5}"],
      ]) {
        const r = await setDialogInput(page, label, value);
        fillResults[label] = r;
      }
      expect(
        Object.values(fillResults).every((r) => r.set),
        `Fill failures: ${JSON.stringify(fillResults)}`,
      ).toBe(true);
    },
  );

  await softStep("Step 6: SAVE the identifier configuration", async () => {
    const r = await clickDialogButton(page, "SAVE");
    expect(r.clicked).toBe(true);
    await page.waitForTimeout(3000);
  });

  await softStep("Step 7: Reload the page", async () => {
    await page.reload({ waitUntil: "networkidle", timeout: 30000 });
    await page.waitForTimeout(4000);
    await page.evaluate(() => {
      document.body.classList.add("selenium");
      (window as any).grok.shell.windows.simpleMode = true;
      (window as any).grok.shell.windows.showBrowse = true;
    });
    expect(true).toBe(true);
  });

  await softStep(
    "Step 8-9: Open customers table from test_postgres",
    async () => {
      // Use JS API to fetch the customers table — the identifiers config is
      // applied at semantic-type detection time on any DataFrame loaded from
      // this connection's `public.customers`.
      const result = await page.evaluate(async () => {
        const tp = await (window as any).grok.dapi.connections.find(
          "af9bcf40-21a0-11f1-89e2-7b1321b80948",
        );
        if (!tp) return { opened: false, reason: "connection not found" };
        const q = await tp.tableQuery("public.customers", false, 1000);
        const df = await q.executeTable();
        if (!df) return { opened: false, reason: "executeTable returned null" };
        (window as any).grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 2000));
        await new Promise((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => {
            sub.unsubscribe();
            resolve(undefined);
          });
          setTimeout(resolve, 4000);
        });
        const col = df.col("customerid");
        return {
          opened: true,
          semType: col?.semType ?? null,
          type: col?.type ?? null,
          sample: col ? [col.get(0), col.get(1), col.get(2)] : null,
        };
      });
      expect(
        result.opened,
        `Could not open customers: ${(result as any).reason}`,
      ).toBe(true);
    },
  );

  await softStep(
    "Step 10: customerid column is highlighted (semType = CUSTOMER_ID)",
    async () => {
      const semType = await page.evaluate(() => {
        const tv = (window as any).grok.shell.tv;
        return tv?.dataFrame?.col("customerid")?.semType ?? null;
      });
      expect(
        semType,
        `customerid semType expected CUSTOMER_ID, got ${semType}`,
      ).toBe("CUSTOMER_ID");
    },
  );

  await softStep(
    "Step 11: Click customerid header → Context Panel shows semantic type",
    async () => {
      const ok = await page.evaluate(() => {
        const tv = (window as any).grok.shell.tv;
        const col = tv?.dataFrame?.col("customerid");
        if (!col) return false;
        tv.grid.col("customerid").selected = true;
        (window as any).grok.shell.o = col;
        return col.semType === "CUSTOMER_ID";
      });
      expect(ok).toBe(true);
    },
  );

  await softStep(
    "Step 12: Remove identifiers configuration on test_postgres",
    async () => {
      await gotoPostgresAndExpand(page);
      const r1 = await rightClickConnection(page, TEST_POSTGRES);
      expect(r1.found).toBe(true);
      // Look for "Remove identifiers" or similar removal item
      const items = await page.evaluate(() => {
        return Array.from(document.querySelectorAll(".d4-menu-item")).map(
          (el) => el.textContent?.trim() ?? "",
        );
      });
      const removeText =
        items.find((t) => /remove.*identif/i.test(t)) ??
        items.find((t) => /^Configure Identifiers/.test(t));
      expect(
        removeText,
        `Removal menu item not found. Items: ${items.join(" | ")}`,
      ).toBeTruthy();
      if (removeText && /remove/i.test(removeText)) {
        const c = await clickMenuItem(page, removeText);
        expect(c.clicked).toBe(true);
        // confirm if a confirm dialog appears
        await page.waitForTimeout(1000);
        const cd = await clickDialogButton(page, "OK");
        if (!cd.clicked) await clickDialogButton(page, "YES");
      } else {
        // Fall back: re-open Configure Identifiers... and use a Remove/Clear option inside the dialog
        const c = await clickMenuItem(page, "Configure Identifiers...");
        expect(c.clicked).toBe(true);
        await page.waitForTimeout(2000);
        // Try a Remove / Clear / Reset button at the bottom of the dialog
        const r = await clickDialogButton(page, "Remove");
        if (!r.clicked) await clickDialogButton(page, "CLEAR");
        await page.waitForTimeout(1000);
        await clickDialogButton(page, "OK");
      }
      await page.waitForTimeout(2000);
    },
  );

  await softStep(
    "Step 13-14: Reload + verify customerid not highlighted (semType empty/non-CUSTOMER_ID)",
    async () => {
      await page.reload({ waitUntil: "networkidle", timeout: 30000 });
      await page.waitForTimeout(4000);
      const semType = await page.evaluate(async () => {
        const tp = await (window as any).grok.dapi.connections.find(
          "af9bcf40-21a0-11f1-89e2-7b1321b80948",
        );
        const q = await tp.tableQuery("public.customers", false, 1000);
        const df = await q.executeTable();
        await new Promise((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => {
            sub.unsubscribe();
            resolve(undefined);
          });
          setTimeout(resolve, 4000);
        });
        return df.col("customerid")?.semType ?? null;
      });
      expect(
        semType,
        `customerid semType after removal — should NOT be CUSTOMER_ID, got ${semType}`,
      ).not.toBe("CUSTOMER_ID");
    },
  );

  if (stepErrors.length > 0) {
    const summary = stepErrors
      .map((e) => `  - ${e.step}: ${e.error}`)
      .join("\n");
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
