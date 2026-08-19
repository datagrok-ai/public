import {Page, expect} from '@playwright/test';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

// Monomer-library selection lives in per-user server settings, so the specs that toggle libraries
// cannot run next to each other: one spec's excluded library breaks the other's detection or
// checkbox state. Playwright serialises only within a file, so take a lock across worker processes.
const LOCK_DIR = path.join(os.tmpdir(), 'dg-bio-monomer-lib.lock');
const STALE_MS = 5 * 60_000;
const WAIT_MS = 4 * 60_000;

export async function acquireMonomerLibLock(): Promise<void> {
  const deadline = Date.now() + WAIT_MS;
  for (;;) {
    try {
      fs.mkdirSync(LOCK_DIR);
      return;
    }
    catch (e) {
      let ageMs = 0;
      try { ageMs = Date.now() - fs.statSync(LOCK_DIR).mtimeMs; }
      catch (_) { continue; }
      // a lock left behind by a killed run must not wedge every later run
      if (ageMs > STALE_MS) {
        fs.rmSync(LOCK_DIR, {recursive: true, force: true});
        continue;
      }
      if (Date.now() > deadline)
        throw new Error(`monomer-library lock held for over ${WAIT_MS / 1000}s by another spec`);
      await new Promise((r) => setTimeout(r, 1000));
    }
  }
}

export function releaseMonomerLibLock(): void {
  fs.rmSync(LOCK_DIR, {recursive: true, force: true});
}

/** Expands a Bio top-menu group. `hover()` jumps straight to the element centre and d4 expands a
 * submenu off the mouse path, so a jump leaves the container at display:none — walk the pointer in
 * instead. Verified on dev: a stepped path expands all six groups within one open menu. */
export async function openBioGroup(page: Page, group: string): Promise<void> {
  const groupLoc = page.locator(`[name="div-Bio---${group}"]`);
  // Clicking the menu bar item toggles, so a second click on an already-open menu closes it.
  if (!await groupLoc.isVisible().catch(() => false))
    await page.locator('[name="div-Bio"]').click();
  await groupLoc.waitFor({state: 'visible', timeout: 15_000});
  const box = (await groupLoc.boundingBox())!;
  await page.mouse.move(box.x + 5, box.y + box.height / 2, {steps: 12});
  await page.mouse.move(box.x + box.width - 15, box.y + box.height / 2, {steps: 12});
  await page.locator(`[name^="div-Bio---${group}---"]`).first()
    .waitFor({state: 'visible', timeout: 15_000});
}

/** Opens a Bio top-menu leaf. Unlike the shared `bio.openBioAnalyze` helper this walks any group,
 * and uses real clicks/hovers — d4 menus ignore synthetic mouseover events. */
export async function openBioMenu(page: Page, group: string, leaf: string): Promise<void> {
  await openBioGroup(page, group);
  const leafLoc = page.locator(`[name="div-Bio---${group}---${leaf}"]`);
  await leafLoc.waitFor({state: 'visible', timeout: 15_000});
  await leafLoc.click();
}

/** Waits for a d4 dialog by its `name` attribute and returns its locator. */
export function dialog(page: Page, name: string) {
  return page.locator(`[name="dialog-${name}"]`);
}

/** Cancels the named dialog and waits for it to detach — a dialog left open intercepts
 * every later click. */
export async function cancelDialog(page: Page, name: string): Promise<void> {
  const dlg = dialog(page, name);
  if (await dlg.count() === 0) return;
  const cancel = dlg.locator('[name="button-CANCEL"]');
  if (await cancel.count() > 0) await cancel.click();
  else await page.keyboard.press('Escape');
  await dlg.waitFor({state: 'detached', timeout: 15_000});
}

/** Sets a d4 choice input inside a dialog, verifying the selection stuck. */
export async function setChoice(page: Page, dlgName: string, input: string, value: string): Promise<void> {
  const sel = page.locator(`[name="dialog-${dlgName}"] [name="input-host-${input}"] select`);
  await sel.selectOption({label: value});
  await expect(sel).toHaveValue(value);
}

/** Types into a d4 text input and blurs — `fill()` updates the DOM but not the Dart value. */
export async function setText(page: Page, dlgName: string, input: string, value: string): Promise<void> {
  const inp = page.locator(`[name="dialog-${dlgName}"] [name="input-host-${input}"] input`);
  await inp.click();
  await inp.press('Control+a');
  await inp.pressSequentially(value);
  await inp.press('Tab');
}

/** Toggles a checkbox input to the requested state (no-op when already there). */
export async function setBool(page: Page, dlgName: string, input: string, value: boolean): Promise<void> {
  const box = page.locator(`[name="dialog-${dlgName}"] [name="input-host-${input}"] input[type="checkbox"]`);
  if (await box.isChecked() !== value) await box.click();
  await expect(box).toBeChecked({checked: value});
}

/** Loads a CSV into a fresh table view, optionally keeping only the first `rows` rows, and
 * waits until Bio has detected the Macromolecule columns and the grid has painted. */
export async function loadBioTable(page: Page, path: string, rows?: number, name?: string): Promise<void> {
  await page.evaluate(async ({p, r, n}) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const src = await grok.dapi.files.readCsv(p);
    const df = r != null && r < src.rowCount ?
      src.clone(DG.BitSet.create(src.rowCount, (i) => i < r)) : src;
    df.name = n ?? df.name;
    grok.shell.addTableView(df);
    const detected = () => df.columns.toList().some((c: any) => c.semType === 'Macromolecule');
    for (let i = 0; i < 200; i++) {
      if (detected() && document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((res) => setTimeout(res, 200));
    }
  }, {p: path, r: rows, n: name});
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.locator('[name="div-Bio"]').waitFor({state: 'visible', timeout: 60_000});
  await page.evaluate(async () => { await (grok as any).functions.call('Bio:getSeqHelper', {}); });
}

/** Closes every viewer the test docked, leaving the grid — keeps specs independent. */
export async function closeExtraViewers(page: Page): Promise<void> {
  await page.evaluate(() => {
    for (const v of Array.from(((grok.shell.tv as any)?.viewers ?? []) as any[]))
      if (v.type !== 'Grid') v.close();
  });
}
