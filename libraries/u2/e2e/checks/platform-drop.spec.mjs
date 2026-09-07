/* P4.5 — the platform's own drag channel: a file or a query dragged out of Browse becomes a data
   source on the tray. The drag is REAL — the source element is built in the page with the very call
   a Browse tree row makes (`ui.makeDraggable`, tree_view.dart:94-98) and driven with the mouse, so
   everything between it and our `mouseup` is the platform's own code: the AppEvents bus, the ghost,
   the dart-side droppable registration. Local mode answers file reads from `api.json`, so a
   dropped file is asserted for a readable ANSWER, never for data — that proof is the stand's. */
import {ok, shot} from '../local.mjs';
import {balloons, chipMenu, chipNames, clearBalloons, components, dragChrome, focusTree, newForm,
  platformDrag, selectChip, setPanelField, spyRuns, surfaceNode, toMode} from '../lib.mjs';

export async function fixture(page) {
  await newForm(page);
}

/** The chip's own Refresh, waited out — the one design-time run a non-query ever gets. */
async function refresh(page, name) {
  await selectChip(page, name);
  await chipMenu(page, name, 'Refresh');
  await page.waitForFunction((name) => [...document.querySelectorAll('.d4-balloon-content, .d4-balloon')]
    .some((b) => b.textContent.startsWith(`${name}:`)), name, {timeout: 8000}).catch(() => {});
}

async function checkPlatformDrop(page) {
  const csv = {path: 'System:AppData/U2demo/orders.csv'};

  // a: the canvas — the wiring, and an answer instead of silence
  await clearBalloons(page);
  await platformDrag(page, csv, page.locator('.u2-designer-surface'));
  const dropped = await components(page);
  const after = await dragChrome(page);
  const said = await balloons(page);
  await shot(page, 'platform-drop-1-file');
  ok('platform-drop/1a/a-file-dropped-on-the-canvas-becomes-an-openfile-source-over-its-own-path',
    JSON.stringify(await chipNames(page)) === JSON.stringify(['orders1']) &&
    dropped.length === 1 && dropped[0].props?.func === 'OpenFile' &&
    dropped[0].props?.params?.fullPath === csv.path &&
    after.body === false && after.ghosts === 0 && said.some((t) => /orders1/.test(t)),
    `${JSON.stringify(dropped)} chrome=${JSON.stringify(after)} said=${JSON.stringify(said)}`);

  await focusTree(page, 'form1');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);
  const undone = await chipNames(page);

  // b: the tray is a drop target too, and says so for the whole drag
  let during = null;
  await platformDrag(page, csv, page.locator('.u2-designer-tray'),
    async () => during = await dragChrome(page));
  await shot(page, 'platform-drop-1-tray');
  ok('platform-drop/1b/the-tray-says-where-the-drop-will-land-and-a-drop-on-it-lands',
    undone.length === 0 && during.body === true && during.tray === true &&
    during.hover === 'none' && (await dragChrome(page)).tray === false &&
    JSON.stringify(await chipNames(page)) === JSON.stringify(['orders1']),
    `undone=${JSON.stringify(undone)} during=${JSON.stringify(during)} ` +
    `chips=${JSON.stringify(await chipNames(page))}`);

  await focusTree(page, 'form1');
  await page.keyboard.press('Control+z');
  await page.waitForTimeout(300);

  // c: a function, dropped on a rendered control — inserted, run by nothing, run by Refresh
  await spyRuns(page);
  await platformDrag(page, {func: 'demoOrders'}, surfaceNode(page, 'form1'));
  const onDrop = await page.evaluate(() => window.__u2Runs);
  const wired = await components(page);
  await clearBalloons(page);
  await refresh(page, 'demoOrders1');
  const empty = await balloons(page);
  await page.evaluate(() => window.__u2Sub?.unsubscribe?.());
  await shot(page, 'platform-drop-1-func');
  ok('platform-drop/1c/a-dropped-function-lands-wired-and-nothing-runs-it',
    wired.length === 1 && wired[0].props?.func === 'demoOrders' &&
    wired[0].props?.params === undefined && !onDrop.includes('demoOrders') &&
    empty.some((t) => /demoOrders1: ready · 0 rows/.test(t)),
    `${JSON.stringify(wired)} runs=${JSON.stringify(onDrop)} said=${JSON.stringify(empty)}`);

  // OR-D2: a drop never opens a dialog to ask for the parameters — the panel's Parameters section
  // is where they are filled in, and the source says what it answered before and after
  await clearBalloons(page);
  await setPanelField(page, 'days', '30');
  await refresh(page, 'demoOrders1');
  const filled = await balloons(page);
  await shot(page, 'platform-drop-1-params');
  ok('platform-drop/1c2/the-parameters-section-is-where-a-dropped-functions-params-are-filled-in',
    filled.some((t) => /demoOrders1: ready · 4 rows/.test(t)) &&
    JSON.stringify((await components(page))[0].props?.params) === JSON.stringify({days: 30}),
    `said=${JSON.stringify(filled)} ${JSON.stringify(await components(page))}`);

  // d: what nothing can open is refused BY NAME — accepting the drag is what makes it visible
  await clearBalloons(page);
  await platformDrag(page, {path: 'System:AppData/U2demo/report.pdf'},
    page.locator('.u2-designer-surface'));
  const refusal = await balloons(page);
  await shot(page, 'platform-drop-1-refused');
  ok('platform-drop/1d/a-file-nothing-opens-as-a-table-is-refused-by-name',
    JSON.stringify(await chipNames(page)) === JSON.stringify(['demoOrders1']) &&
    refusal.some((t) => /report\.pdf/.test(t) && /pdf/.test(t)),
    `chips=${JSON.stringify(await chipNames(page))} said=${JSON.stringify(refusal)}`);

  // e: Run mode has no tray to point at, so the drop is refused before it starts
  await toMode(page, 'Run');
  let running = null;
  await platformDrag(page, csv, page.locator('.u2-designer-surface'),
    async () => running = await dragChrome(page));
  await toMode(page, 'Design');
  await page.waitForTimeout(300);
  ok('platform-drop/1e/a-drop-in-run-mode-adds-nothing-and-highlights-nothing',
    running.tray === false &&
    JSON.stringify(await chipNames(page)) === JSON.stringify(['demoOrders1']),
    `during=${JSON.stringify(running)} chips=${JSON.stringify(await chipNames(page))}`);
}

export const checks = [
  {id: 'platform-drop/1 platform drop: a dragged file and a dragged function become sources, refusals by name',
    run: checkPlatformDrop},
];
