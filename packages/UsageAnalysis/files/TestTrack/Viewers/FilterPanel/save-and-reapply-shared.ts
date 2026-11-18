import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';
import {HEADER_COUNTER} from '../../helpers/filter-panel';
import {isLocalBootNoise, softStep} from '../../spec-login';

declare const grok: any;

const AMBIENT_CONSOLE_ERROR = 'Permissions policy violation: compute-pressure';

export async function readCounterTooltipSummary(page: Page): Promise<Record<string, string> | null> {
  return page.evaluate(async (sel) => {
    const counterEl = document.querySelector(sel);
    if (!counterEl) return null;
    const r = counterEl.getBoundingClientRect();
    const at = {clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, bubbles: true};
    counterEl.dispatchEvent(new MouseEvent('mouseenter', at));
    counterEl.dispatchEvent(new MouseEvent('mousemove', at));
    counterEl.dispatchEvent(new MouseEvent('mouseover', at));
    const table = await (window as any).__poll(() => document.querySelector('.d4-tooltip table'),
      (t: Element | null) => t !== null, 2000, 50) as Element | null;
    if (!table) return null;
    const summary: Record<string, string> = {};
    for (const row of Array.from(table.querySelectorAll('tr'))) {
      const cells = Array.from(row.querySelectorAll('td,th'));
      if (cells.length < 2) continue;
      const caption = cells[0].textContent?.trim() ?? '';
      if (caption.length > 0) summary[caption] = cells[1].textContent?.trim() ?? '';
    }
    return Object.keys(summary).length > 0 ? summary : null;
  }, HEADER_COUNTER);
}

export async function readSaveOrApplyLeaves(page: Page): Promise<string[]> {
  await page.evaluate(() => {
    const root = grok.shell.tv.root as HTMLElement;
    const filtersViewer = root.querySelector('[name="viewer-Filters"]');
    const titlebar = filtersViewer?.closest('.panel-base')?.querySelector('.panel-titlebar');
    const hamburger = titlebar?.querySelector('[name="icon-font-icon-menu"]') as HTMLElement | null;
    if (!hamburger) {
      throw new Error('no [name="icon-font-icon-menu"] in the Filter Panel title bar of view ' +
        `"${grok.shell.tv.name}"`);
    }
    hamburger.click();
  });
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 15_000});
  return page.evaluate(async () => {
    const w = window as any;
    const popup = () => Array.from(document.querySelectorAll('.d4-menu-popup')).pop() ?? null;
    const labelsOf = (root: Element | null) =>
      Array.from(root?.querySelectorAll('.d4-menu-item-label') ?? [])
        .map((l) => (l.textContent ?? '').trim());
    const findGroup = () => Array.from(popup()?.children ?? [])
      .filter((c) => c.classList.contains('d4-menu-item'))
      .find((it) => it.querySelector(':scope > .d4-menu-item-label')?.textContent?.trim() ===
        'Save or Apply') ?? null;
    const group = await w.__poll(findGroup, (g: Element | null) => g !== null, 5000, 50) as Element | null;
    if (!group) {
      throw new Error('panel menu: group "Save or Apply" not found in the panel popup; visible: ' +
        labelsOf(popup()).join(' | '));
    }
    const rect = group.getBoundingClientRect();
    const at = {clientX: rect.x + rect.width / 2, clientY: rect.y + rect.height / 2, bubbles: true};
    for (const type of ['mouseover', 'mouseenter', 'mousemove'])
      group.dispatchEvent(new MouseEvent(type, at));
    const leaves = () => Array.from(group!.querySelectorAll('.d4-menu-item'))
      .map((c) => c.querySelector(':scope > .d4-menu-item-label')?.textContent?.trim() ?? '')
      .filter((s) => s.length > 0);
    return w.__poll(leaves, (out: string[]) => out.length > 0, 5000, 50);
  });
}

async function visibleMenuPopups(page: Page): Promise<number> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.d4-menu-popup')).filter((p) => {
      const cs = window.getComputedStyle(p as HTMLElement);
      const r = (p as HTMLElement).getBoundingClientRect();
      return cs.display !== 'none' && cs.visibility !== 'hidden' && r.width > 0 && r.height > 0;
    }).length);
}

export async function dismissPanelMenu(page: Page): Promise<void> {
  for (let attempt = 0; attempt < 3; attempt++) {
    await page.keyboard.press('Escape');
    if (await v.pollValue(() => visibleMenuPopups(page), (n) => n === 0, 3000, 40) === 0) return;
  }
  expect(await visibleMenuPopups(page),
    'the panel menu popup stayed on screen and would intercept the next gesture').toBe(0);
}

// The caller's own saved filter presets of that name: the ids of the SavedFilter entities on the server.
const myPresetIds = (page: Page, name: string): Promise<string[]> => page.evaluate(async (n) => {
  const headers = {Authorization: String(grok.dapi.token)};
  const found: any[] = await (await fetch(`${grok.dapi.root}/entities?text=${encodeURIComponent(n)}`, {headers})).json();
  return found.filter((e) => e['#type'] === 'SavedFilter' && e.friendlyName === n && e.author?.id === grok.shell.user.id)
    .map((e) => e.id as string);
}, name);

// Step 4 of the scenario: the hamburger's Save or Apply > Save... dialog, checked against the
// preset the product saves on the server and against the console.
export async function saveStateViaMenu(page: Page, probeName: string): Promise<void> {
  const consoleErrors: string[] = [];
  const onConsole = (msg: import('@playwright/test').ConsoleMessage) => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  };
  page.on('console', onConsole);
  const errorSet = () => new Set(consoleErrors.filter((t) => !t.includes(AMBIENT_CONSOLE_ERROR) && !isLocalBootNoise(t)));
  const errorsBefore = errorSet();

  await v.drivePanelMenuLeaf(page, 'Filters', 'Save or Apply', 'Save...');

  const dlg = page.locator('.d4-dialog[name="dialog-Save-filter-preset"]');
  await dlg.waitFor({timeout: 10_000});
  const nameInput = dlg.locator('input[name="input-Name"]');
  // Modal.editValue defers focus()+select() to a Timer.run (d4 widgets/dialog/dialog.dart:321).
  // Under two workers that timer lands in the middle of a per-character type, and the next
  // character replaces the whole selection — the dialog then saves under a truncated name, which
  // is why the probe was never saved. Wait the timer out, then set the name in one
  // assignment so a later select() has nothing left to swallow.
  await page.evaluate(() => (window as any).__poll(() => {
    const el = document.querySelector('.d4-dialog[name="dialog-Save-filter-preset"] input[name="input-Name"]') as
      HTMLInputElement | null;
    return !!el && document.activeElement === el && el.selectionStart === 0 &&
      el.selectionEnd === (el.value ?? '').length;
  }, (ready: boolean) => ready, 1200, 25));
  await nameInput.fill(probeName);
  expect(await v.pollValue(() => nameInput.inputValue(), (val) => val === probeName, 2000, 50),
    'the Save filter preset dialog\'s Name box does not hold the probe name, so OK would save the state ' +
    'under some other name').toBe(probeName);
  await dlg.locator('[name="button-OK"]').click();
  await expect.poll(async () => dlg.count(), {
    timeout: 10_000,
    intervals: [30, 60, 120, 250, 500, 1000],
    message: 'the Save filter preset dialog did not close after OK',
  }).toBe(0);

  // the preset is saved a beat after the dialog closes, so a one-shot read races it under load
  const stored = await v.pollValue(() => myPresetIds(page, probeName), (ids) => ids.length > 0, 10000, 250);
  expect(stored.length, 'no saved filter preset with the probe name is on the server — the state was not saved')
    .toBeGreaterThan(0);

  page.off('console', onConsole);
  const newErrors = [...errorSet()].filter((t) => !errorsBefore.has(t));
  expect(newErrors, `the save produced new console error texts: ${newErrors.join(' | ')}`).toEqual([]);
}

export async function removeProbeState(page: Page, probeName: string): Promise<void> {
  for (const id of await myPresetIds(page, probeName))
    await page.evaluate(async (i) => {
      await fetch(`${grok.dapi.root}/domains/filters/${i}`, {method: 'DELETE', headers: {Authorization: String(grok.dapi.token)}});
    }, id);
  await softStep('Teardown probe preset removed from the server', async () => {
    expect(await myPresetIds(page, probeName), 'the probe preset is still on the server and would bleed into later runs')
      .toEqual([]);
  });
}
