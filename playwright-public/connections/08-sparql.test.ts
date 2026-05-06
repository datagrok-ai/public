import { test, expect } from '@playwright/test';
import {
  AUTH_STATE,
  SPARQL_ENDPOINT,
  applyAutomationSetup,
  clickConnectionOk,
  clickMenuItemExact,
  connectionNodeName,
  deleteConnectionByFriendlyName,
  deleteTreeNodeViaContext,
  expandDbProvider,
  fillConnectionField,
  findConnectionByFriendlyName,
  goHome,
} from './helpers';

// Manual scenario `sparql.md` (order 8).
//
// 1. Browse > Databases. Press three dots ("...") to list data sources without connections.
// 2. Right-click Sparql → Add new connection. (UI shortcut)
// 3. Name = test_sparql
// 4. Endpoint, Requires Server = true, Prefixes empty
// 5. Click Test
// 6. Click OK
// 7. Right-click test_sparql → Delete, YES
//
// SCOPE NOTE: on dev the `Show more...` footer under Databases doesn't reveal
// the Sparql provider node directly via a tree click — providers without saved
// connections aren't surfaced as right-clickable nodes. Instead we use the
// equivalent UI affordance: the global **New Connection** button surfaces a
// dialog with a `Data Source` dropdown that lists *all* providers (including
// Sparql). This matches the existing `sparql-spec.ts` template which drives
// the same path. The Sparql dialog also has no TEST button on dev (per the
// existing spec), so step 5 is omitted from the autotest — see sparql-ui.md
// for a manual TEST-button pass when an interactive Sparql endpoint is up.

const PROVIDER = 'Sparql';
const CONNECTION = 'test_sparql';

test.describe.serial('Connections / SPARQL', () => {
  test.beforeAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, CONNECTION);
    await ctx.close();
  });

  test.afterAll(async ({ browser }) => {
    const ctx = await browser.newContext({ storageState: AUTH_STATE });
    const page = await ctx.newPage();
    await goHome(page);
    await deleteConnectionByFriendlyName(page, CONNECTION);
    await ctx.close();
  });

  test('1. Add Sparql via New Connection dialog; OK persists test_sparql', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    // The New Connection button lives on the Databases browse view (`/db?browse=db`),
    // not on home. Navigate explicitly so the button is rendered.
    await page.goto('/db?browse=db', { waitUntil: 'domcontentloaded' });
    // The Connections view paints its toolbar (where the NEW CONNECTION button
    // lives) lazily — poll for the button to attach instead of relying on a
    // fixed wait. On dev the paint can take 5–15s after a fresh navigation,
    // especially after several heavy preceding tests have touched the tree.
    const clicked = await page
      .waitForFunction(() => {
        const btn = (Array.from(document.querySelectorAll('button, .ui-btn')) as HTMLElement[])
          .find((b) => (b.textContent ?? '').trim().toUpperCase().startsWith('NEW CONNECTION')
            && b.offsetParent !== null);
        if (!btn) return false;
        btn.click();
        return true;
      }, undefined, { timeout: 30_000 })
      .then(() => true)
      .catch(() => false);
    expect(clicked, 'New Connection button visible').toBe(true);

    // Wait for the named dialog and the Data Source select to actually be in the DOM —
    // the platform paints the form lazily after the global click.
    await page.locator('[name="dialog-Add-new-connection"]').waitFor({ timeout: 10_000 });
    await page.waitForSelector(
      '[name="dialog-Add-new-connection"] [name="input-host-Data-Source"] select',
      { state: 'attached', timeout: 10_000 },
    );

    // Change `Data Source` to Sparql via Playwright's native selectOption — that
    // dispatches the full `input`/`change` pair the way real keyboard/mouse input
    // does, which Dart input listeners on the choice widget rely on. JS `.value`
    // setter alone wasn't enough on dev (the Sparql-specific form didn't render).
    await page.locator('[name="dialog-Add-new-connection"] [name="input-host-Data-Source"] select')
      .selectOption(PROVIDER);

    // Wait for the Sparql-specific input layout to render — change-event triggers
    // a full form re-render which is slow on dev (Dart-side).
    await page.waitForSelector(
      '[name="dialog-Add-new-connection"] [name="input-host-Endpoint"]',
      { state: 'visible', timeout: 30_000 },
    );

    await fillConnectionField(page, 'Name', CONNECTION);
    await fillConnectionField(page, 'Endpoint', SPARQL_ENDPOINT);

    // "Requires Server" boolean checkbox — toggle if not already on.
    await page.evaluate(() => {
      const cb = document.querySelector(
        '.d4-dialog [name="input-host-Requires-Server"] input[type="checkbox"]',
      ) as HTMLInputElement | null;
      if (cb && !cb.checked) cb.click();
    });

    // Sparql dialog has no TEST button on dev — go straight to OK.
    await clickConnectionOk(page);

    const saved = await findConnectionByFriendlyName(page, CONNECTION);
    expect(saved, `connection "${CONNECTION}" should exist after OK`).not.toBeNull();
    expect(saved!.dataSource).toBe(PROVIDER);
  });

  test('2. Delete test_sparql via context menu', async ({ page }) => {
    await goHome(page);
    await applyAutomationSetup(page);
    await expandDbProvider(page, PROVIDER);

    await deleteTreeNodeViaContext(page, connectionNodeName(PROVIDER, CONNECTION));

    await expect.poll(async () => (await findConnectionByFriendlyName(page, CONNECTION)) === null,
      { timeout: 15_000 }).toBe(true);
  });
});
