import {Page} from '@playwright/test';

declare const grok: any;

// Section-local speed helpers. Same actuations and same guards as helpers/chem.ts and
// spec-login.ts, with the fixed delays replaced by polls on the thing that actually changes.

// helpers/chem.ts openChemMenuItem sleeps a flat 600 ms between opening the Chem menubar root and
// looking for the leaf. The leaf usually lands in the document within one frame, so the sleep is
// pure wall clock; poll for it instead, keeping the same cap and the same owner check.
export async function openChemMenuItemFast(
  page: Page, label: string, opts?: {delayMs?: number},
): Promise<void> {
  const capMs = opts?.delayMs ?? 600;
  await page.evaluate(async ({label, capMs}) => {
    const chemMenu = document.querySelector('[name="div-Chem"]') as HTMLElement;
    // Menu labels from a previously opened menu stay in the document. Polling for the label alone
    // therefore returns instantly with a STALE node, and clicking it actuates nothing — which is
    // what the flat 600 ms sleep was hiding by outlasting the old menu's teardown. Only a node
    // that was not already there counts as this menu's leaf. The root also TOGGLES, so a menu left
    // open means the first click closes it; retry once.
    const stale = new Set(Array.from(document.querySelectorAll('.d4-menu-item-label')));
    const findFresh = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
      .find((m) => !stale.has(m) && m.textContent!.trim() === label) as HTMLElement | undefined;
    let item: HTMLElement | undefined;
    for (let attempt = 0; attempt < 2 && !item; attempt++) {
      chemMenu.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      const deadline = Date.now() + capMs;
      item = findFresh();
      while (!item && Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 25));
        item = findFresh();
      }
    }
    // Last resort: the same node the flat-sleep version would have taken.
    if (!item)
      item = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .find((m) => m.textContent!.trim() === label) as HTMLElement | undefined;
    if (!item) throw new Error(`openChemMenuItemFast("${label}"): leaf never appeared within ${capMs} ms`);
    const menuItem = item.closest('.d4-menu-item') as HTMLElement;
    // One click on [name="div-Chem"] puts every menubar root's labels in the document, so an exact
    // text match can land on another feature's item; keep helpers/chem.ts' owner guard verbatim.
    const named = item.closest('[name^="div-"]');
    const owner = named ? named.getAttribute('name') : null;
    if (owner !== null && !owner.startsWith('div-Chem'))
      throw new Error(`openChemMenuItemFast("${label}") matched ${owner}, which is not a Chem menu item`);
    menuItem.dispatchEvent(new MouseEvent('click', {bubbles: true}));
  }, {label, capMs});
}

// spec-login.ts waitForMolecule only ever resolves on semType 'Molecule'. A table whose molecular
// column carries another semantic type (ChemicalMixture) can never satisfy it, so the call burns
// its whole 45 s cap and the caller's `.catch(() => {})` hides it. This is the same wait, keyed on
// the semantic type the table under test actually declares.
export async function waitForSemType(page: Page, semType: string, timeoutMs = 45_000): Promise<void> {
  await page.evaluate(({st, timeout}) => new Promise<void>((resolve, reject) => {
    const g = (window as any).grok;
    const typed = () => [g?.shell?.t, (window as any).__df]
      .filter(Boolean)
      .some((t: any) => t.columns.toList().some((c: any) => c.semType === st));

    if (typed())
      return resolve();

    const sub = g.events.onEvent('ddt-semantic-type-detected').subscribe(() => { if (typed()) done(); });
    const poll = setInterval(() => { if (typed()) done(); }, 100);
    const timer = setTimeout(
      () => done(new Error(`waitForSemType: no ${st} column detected within ${timeout} ms`)), timeout);

    function done(err?: Error) {
      clearTimeout(timer);
      clearInterval(poll);
      sub.unsubscribe();
      err ? reject(err) : resolve();
    }
  }), {st: semType, timeout: timeoutMs});
}

// The context panel rebuilds its accordion asynchronously after grok.shell.o is set. The specs
// stood a flat sleep in for that; poll the pane set until it stops changing instead, capped at the
// sleep replaced. A panel that never settles still burns the cap and fails the caller's assertion.
export async function settleContextPanes(page: Page, capMs: number): Promise<void> {
  await page.evaluate(async (cap: number) => {
    const read = () => Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
      .map((h) => (h.textContent ?? '').trim()).join('|');
    const deadline = Date.now() + cap;
    let last = '';
    let stable = 0;
    while (Date.now() < deadline) {
      const now = read();
      if (now.length > 0 && now === last) { if (++stable >= 3) return; }
      else stable = 0;
      last = now;
      await new Promise((r) => setTimeout(r, 100));
    }
  }, capMs);
}

// The grid keeps repainting after its canvas attaches, which the specs stood a flat sleep in for.
// Hold the canvas's own pixels steady instead, capped at the sleep replaced.
export async function settleGridPaint(page: Page, capMs: number): Promise<void> {
  await page.evaluate(async (cap: number) => {
    const deadline = Date.now() + cap;
    let last = -1;
    let stable = 0;
    while (Date.now() < deadline) {
      const cv = document.querySelector('[name="viewer-Grid"] canvas') as HTMLCanvasElement | null;
      let h = -1;
      if (cv && cv.width > 0 && cv.height > 0) {
        const d = cv.getContext('2d')!.getImageData(0, 0, cv.width, Math.min(cv.height, 160)).data;
        h = 0;
        for (let i = 0; i < d.length; i += 997) h = (h * 31 + d[i]) | 0;
      }
      if (h !== -1 && h === last) { if (++stable >= 3) return; }
      else stable = 0;
      last = h;
      await new Promise((r) => setTimeout(r, 100));
    }
  }, capMs);
}

// A bounded wait whose expiry is not itself a failure: it stands in for a flat sleep of the same
// length, so reaching the cap must leave the caller's own assertion to decide.
export async function waitQuiet(p: Promise<unknown>): Promise<void> {
  try { await p; }
  catch (_) { }
}

// The Chem menubar root landing in the DOM — what the flat post-login settles stood in for.
export async function waitForChemMenuRoot(page: Page, capMs = 3000): Promise<void> {
  await waitQuiet(page.waitForFunction(() => !!document.querySelector('[name="div-Chem"]'),
    null, {timeout: capMs, polling: 100}));
}
