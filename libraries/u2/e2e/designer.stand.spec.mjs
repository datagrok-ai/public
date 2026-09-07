/* Lane 2 — what only a stand can answer: that the package as published (not as staged) registers
   the app, that the platform's own routes reach it, and that it is filed where Browse expects it.
   Everything about the designer's behaviour is Lane 1's job — do not grow this file into a copy of
   checks/. Run with `npm run e2e:stand` on a stand nobody else is using; publish first
   (`npm run debug-u2demo-<host>` in packages/U2Demo). */
import {ARTIFACTS, note, ok, report, results, shot} from './local.mjs';
import {balloons, bindThroughPicker, clearBalloons, components, dropControl, fieldValue,
  newForm, openSpec, panel, platformDrag, selectChip, specErrors, toMode, viewerAt, waitStatus}
  from './lib.mjs';
import {STAND_URL, openAppByRoute, openAppFromBrowse, openStand} from './stand.mjs';

const APP = {package: 'U2demo', route: 'U2Designer', friendlyName: 'U2 Designer', browsePath: 'Dev'};

const designerState = (page) => page.evaluate(() => ({
  designer: !!document.querySelector('.u2-designer'),
  view: grok.shell.v?.name,
  toolbox: !!document.querySelector('.u2-designer-toolbox .u2-tree'),
  named: [...document.querySelectorAll('.u2-designer-surface [data-u2-name]')].map((e) => e.dataset.u2Name),
  status: [...document.querySelectorAll('div')]
    .find((x) => x.children.length === 0 && /\d+ node/.test(x.textContent))?.textContent,
}));

/** The published package registers the app as a function, with the metadata the client routes on. */
async function checkRegistration(page) {
  const found = await page.evaluate((pkg) => DG.Func.find({package: pkg, name: 'u2DesignerApp'})
    .map((f) => ({name: f.name, friendlyName: f.friendlyName, tags: f.tags,
      browsePath: f.options?.browsePath})), APP.package);
  ok('S1/app-func-published', found.length === 1 && found[0].browsePath === APP.browsePath,
    JSON.stringify(found));
}

/** The URL route only resolves through the published package — and it is what a shared link is. */
async function checkRoute(page) {
  const opened = await openAppByRoute(page, APP.package, APP.route, '.u2-designer');
  const state = await designerState(page);
  await shot(page, 'stand-1-route');
  ok('S2/url-route-opens-designer', opened && state.designer && /\d+ nodes/.test(state.status ?? ''),
    JSON.stringify(state));
  ok('S3/published-build-is-current', state.named.includes('saveButton'),
    `data-u2-name on the published build: ${JSON.stringify(state.named)}`);
}

/** `meta.browsePath: Dev` — the app has to be reachable without knowing its URL. */
async function checkBrowse(page) {
  await openAppFromBrowse(page, APP.browsePath, APP.friendlyName);
  const state = await designerState(page);
  await shot(page, 'stand-2-browse');
  ok('S4/browse-dev-card-opens-designer', state.designer && state.toolbox, JSON.stringify(state));
}

/* The two data checks canned responses cannot answer (Q10): that a source lists what THIS server
   holds, and that the design-time row cap survives a real query execution. */

const USERS_SPEC = JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-entity-source', name: 'people', props: {entity: 'users', pageSize: 10}}],
  root: {tag: 'u2-form', name: 'usersForm', children: [
    {tag: 'u2-choice-input', name: 'userPick', props: {label: 'User'}, bind: {items: '$.people.names'}},
    {tag: 'u2-text-input', name: 'stateField', props: {label: 'State'}, bind: {value: '$.people.state'}},
  ]},
}, null, 2);

/** No `designData`: a query defaults to `live` (DD9), which is the half of the promise the cap
 * completes — the source runs for real at DESIGN time, and it runs capped. */
const cappedSpec = (func) => JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-func-source', name: 'capped', props: {func, debounce: 50}}],
  root: {tag: 'u2-form', name: 'capForm', children: [
    {tag: 'u2-text-input', name: 'rowsField', props: {label: 'Rows'}, bind: {value: '$.capped.rowCount'}},
    {tag: 'u2-text-input', name: 'stateField', props: {label: 'State'}, bind: {value: '$.capped.state'}},
  ]},
}, null, 2);

/** `u2-entity-source` over the real `grok.dapi.users` — the same spec local mode answers from
 * `api.json` fixtures, here against whoever this server actually has. */
async function checkEntitySource(page) {
  await openSpec(page, USERS_SPEC);
  await waitStatus(page, 'usersForm');
  await page.waitForTimeout(3000);
  const shown = await page.evaluate(() =>
    [...document.querySelectorAll('.u2-designer-surface [data-u2-name="userPick"] option')]
      .map((o) => o.textContent).filter((t) => t !== ''));
  const real = await page.evaluate(async () =>
    (await grok.dapi.users.list({pageSize: 100})).map((u) => u.friendlyName ?? u.name));
  // 'done' is the LAST page loaded; a server with more users than one page settles on 'idle'
  const state = await fieldValue(page, 'stateField');
  await shot(page, 'stand-3-users');
  ok('S5/entity-source-lists-this-servers-users',
    shown.length > 0 && shown.every((name) => real.includes(name)) &&
    (state === 'idle' || state === 'done'),
    `state="${state}" shown=${JSON.stringify(shown)} of ${real.length} on the server`);
}

/** A visual query this stand can cap: a `TableQuery` over the Northwind `orders` table (830 rows),
 * saved under `name` and run once for the truth about its row count. The caller removes it. */
const makeProbeQuery = (page, name) => page.evaluate(async (name) => {
  const conns = await grok.dapi.connections.list();
  const connection = conns.find((c) => c.nqName === 'Admin:NW');
  if (!connection)
    return {error: 'no Admin:NW connection on this stand'};
  const query = DG.TableQuery.create(connection);
  query.table = 'orders';
  query.name = name;
  const saved = await grok.dapi.queries.save(query);
  const found = DG.Func.find({name});
  return {id: saved.id, name: saved.name, found: found.length,
    rows: (await found[0].prepare({}).call()).outputs.result.rowCount};
}, name);

const dropProbeQuery = (page, id) =>
  page.evaluate(async (id) => grok.dapi.queries.delete(await grok.dapi.queries.find(id)), id);

/** The design-time 100-row cap (DD9/OR-1). It rides `TableQuery.limit`, and this stand carries no
 * visual query, so the check makes one over a table with 830 rows and removes it again. */
async function checkDesignRowCap(page) {
  const made = await makeProbeQuery(page, 'u2CapProbe');
  if (made.error !== undefined || made.found !== 1) {
    note('S6/row-cap', `not provable here: ${made.error ?? `find answered ${made.found} funcs`}`);
    return;
  }
  try {
    await openSpec(page, cappedSpec(made.name));
    await waitStatus(page, 'capForm');
    await page.waitForTimeout(6000);
    const rows = await fieldValue(page, 'rowsField');
    const state = await fieldValue(page, 'stateField');
    await shot(page, 'stand-4-row-cap');
    ok('S6/the-design-time-run-is-capped-at-100-rows',
      made.rows > 100 && state === 'ready' && Number(rows) === 100 && await specErrors(page) === 0,
      `${made.name} holds ${made.rows} rows; the designer's live run shows ${rows} (state ${state})`);
  } finally {
    await dropProbeQuery(page, made.id);
  }
}

const viewerCapSpec = (func) => JSON.stringify({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-func-source', name: 'capped', props: {func, designData: 'live', debounce: 50}}],
  root: {tag: 'u2-panel', name: 'capPanel', children: [
    {tag: 'u2-viewer-grid', name: 'grid', bind: {table: '$.capped'}},
  ]},
}, null, 2);

/** The cap reaches a viewer (viewers/plan.md WO-V4): a `u2-viewer-grid` over a live query source
 * shows the capped frame at design time and every row in Run mode. */
async function checkViewerRowCap(page) {
  const made = await makeProbeQuery(page, 'u2ViewerCapProbe');
  if (made.error !== undefined || made.found !== 1) {
    note('S6b/viewer-row-cap', `not provable here: ${made.error ?? `find answered ${made.found} funcs`}`);
    return;
  }
  try {
    await openSpec(page, viewerCapSpec(made.name));
    await waitStatus(page, 'capPanel');
    await page.waitForTimeout(6000);
    const design = await viewerAt(page, 'grid');
    await shot(page, 'stand-4b-viewer-cap-design');
    await toMode(page, 'Run');
    await page.waitForTimeout(6000);
    const run = await viewerAt(page, 'grid');
    await shot(page, 'stand-4b-viewer-cap-run');
    await toMode(page, 'Design');
    ok('S6b/the-design-time-cap-reaches-a-grid-viewer-and-run-mode-shows-every-row',
      made.rows > 100 && design?.type === 'Grid' && design.rows > 0 && design.rows <= 100 &&
      run?.rows === made.rows && await specErrors(page) === 0,
      `${made.name} holds ${made.rows} rows; design grid=${JSON.stringify(design)} run grid=${JSON.stringify(run)}`);
  } finally {
    await dropProbeQuery(page, made.id);
  }
}

/* P4.5 — what a drop is FOR. Lane 1 drives the same real platform drag, but local mode answers file
   reads and query runs from `api.json`, so only a stand can say that the data behind a dropped file
   is this server's data and that a dropped query runs capped. */

/** The first tabular file this stand actually holds, and the truth about it read straight off the
 * server through the very function the drop wires (`OpenFile`) — rows, and a column with a value. */
const pickFile = (page) => page.evaluate(async (dirs) => {
  for (const dir of dirs) {
    const list = await grok.dapi.files.list(dir).catch(() => []);
    const csv = list.filter((f) => !f.isDirectory && f.extension === 'csv');
    // a named favourite where the stand has one: the bind picker's tree is virtualized, so the
    // column the check binds has to be near the top of a short-ish column list
    const file = csv.find((f) => /\/(demog|cars|orders)\.csv$/.test(f.fullPath)) ?? csv[0];
    if (file === undefined)
      continue;
    const df = await grok.functions.call('OpenFile', {fullPath: file.fullPath});
    const col = df.columns.names().find((name) => String(df.get(name, 0) ?? '') !== '');
    return {path: file.fullPath, rows: df.rowCount, col, value: String(df.get(col, 0))};
  }
  return {error: `no csv under ${dirs.join(', ')}`};
}, ['System:DemoFiles/', 'System:AppData/ApiSamples/', 'System:AppData/U2demo/']);

/** S7 — a file dragged out of the platform becomes a source over the platform's own reader, and the
 * one run the drop makes (DD9/Q8) reports THIS server's row count, not a fixture's. */
async function checkDroppedFile(page) {
  const file = await pickFile(page);
  if (file.error !== undefined) {
    note('S7/dropped-file', `not provable here: ${file.error}`);
    return;
  }
  await newForm(page);
  await clearBalloons(page);
  await platformDrag(page, {dapiPath: file.path}, page.locator('.u2-designer-surface'));
  const said = await balloons(page);
  const dropped = await components(page);
  const name = dropped[0]?.name;
  if (name !== undefined)
    await selectChip(page, name);
  const status = (await panel(page)).Status ?? {};
  await shot(page, 'stand-5-dropped-file');
  ok('S7a/a-dropped-file-reads-this-servers-data',
    dropped.length === 1 && dropped[0].props?.func === 'OpenFile' &&
    dropped[0].props?.params?.fullPath === file.path && file.rows > 0 &&
    said.some((t) => t === `${name}: ready · ${file.rows} rows`) &&
    status.State === 'ready' && Number(status.Rows) === file.rows,
    `${file.path} holds ${file.rows} rows; ${JSON.stringify(dropped)} ` +
    `status=${JSON.stringify(status)} said=${JSON.stringify(said)}`);

  await dropControl(page, 'u2-text-input', 'text-input', 'form1');
  await dropControl(page, 'u2-number-input', 'number-input', 'form1');
  const walked = await bindThroughPicker(page, 'textInput1', `${file.col} :`, name,
    'currentRow :');
  await bindThroughPicker(page, 'numberInput1', 'currentRowIdx', name);
  await toMode(page, 'Run');
  await page.waitForTimeout(1500);
  const row = page.locator('.u2-designer-surface [data-u2-name="numberInput1"] input').first();
  await row.fill('0');
  await row.press('Enter');
  await page.waitForTimeout(1200);
  const shown = await fieldValue(page, 'textInput1');
  await shot(page, 'stand-6-dropped-file-bound');
  await toMode(page, 'Design');
  ok('S7b/a-field-bound-to-the-dropped-files-own-column-shows-its-own-value',
    walked.found !== undefined && shown === file.value && await specErrors(page) === 0,
    `${name}.${file.col} row 0 is "${file.value}" on the server, the form shows "${shown}" ` +
    `(picker offered ${JSON.stringify(walked.labels?.slice(-6) ?? [])})`);
}

/** S8 — a real query. A `TableQuery` is the only kind the platform lets us cap, so the check makes
 * one over a table with more rows than the cap and removes it again; the qualified-name rule and the
 * DD9 silence are asserted against whatever ambiguous hand-written query this stand already has. */
async function checkDroppedQuery(page) {
  const made = await makeProbeQuery(page, 'u2DropProbe');
  if (made.error !== undefined || made.found !== 1) {
    note('S8/dropped-query', `not provable here: ${made.error ?? `find answered ${made.found} funcs`}`);
    return;
  }
  try {
    await newForm(page);
    await spyRuns(page);
    await platformDrag(page, {func: made.name}, page.locator('.u2-designer-surface'));
    await page.waitForTimeout(4000);
    const dropped = await components(page);
    const ran = await page.evaluate(() => window.__u2Runs);
    await selectChip(page, dropped[0]?.name ?? made.name);
    const status = (await panel(page)).Status ?? {};
    await shot(page, 'stand-7-dropped-query');
    ok('S8a/a-dropped-table-query-runs-itself-at-design-time-and-runs-capped',
      dropped.length === 1 && dropped[0].props?.func === made.name &&
      dropped[0].props?.params === undefined && ran.includes(made.name) &&
      made.rows > 100 && status.State === 'ready' && Number(status.Rows) === 100,
      `${made.name} holds ${made.rows} rows; ${JSON.stringify(dropped)} ` +
      `status=${JSON.stringify(status)} ran=${JSON.stringify(ran)}`);
  } finally {
    await dropProbeQuery(page, made.id);
  }

  // the WO-16 naming rule and the DD9 gate, on a query this stand already has: hand-written SQL
  // whose bare name would name more than one function
  const ambiguous = await page.evaluate(async () => {
    const queries = await grok.dapi.queries.list({pageSize: 400});
    for (const query of queries) {
      const found = DG.Func.find({name: query.name});
      if (found.length > 1 && !(found[0] instanceof DG.TableQuery))
        return {name: found[0].name, nq: found[0].nqName, matches: found.length};
    }
    return {error: 'no ambiguously-named query on this stand'};
  });
  if (ambiguous.error !== undefined) {
    note('S8b/qualified-name', `not provable here: ${ambiguous.error}`);
    return;
  }
  await newForm(page);
  await spyRuns(page);
  await platformDrag(page, {func: ambiguous.name}, page.locator('.u2-designer-surface'));
  await page.waitForTimeout(3000);
  const dropped = await components(page);
  const ran = await page.evaluate(() => window.__u2Runs);
  await selectChip(page, dropped[0]?.name ?? '');
  const status = (await panel(page)).Status ?? {};
  await shot(page, 'stand-8-dropped-ambiguous-query');
  ok('S8b/an-ambiguous-name-is-qualified-and-hand-written-sql-is-not-run-by-a-drop',
    dropped.length === 1 && dropped[0].props?.func === ambiguous.nq &&
    !ran.includes(ambiguous.name) && (status.State ?? '').startsWith('idle'),
    `${ambiguous.name} names ${ambiguous.matches} functions; the spec says ` +
    `"${dropped[0]?.props?.func}" status=${JSON.stringify(status)} ran=${JSON.stringify(ran)}`);
}

/** Every function the platform runs from here on — the differential DD9 gate, not UI inspection. */
const spyRuns = (page) => page.evaluate(() => {
  window.__u2Sub?.unsubscribe?.();
  window.__u2Runs = [];
  window.__u2Sub = grok.functions.onBeforeRunAction
    .subscribe((call) => window.__u2Runs.push(call?.func?.name ?? '?'));
});

const started = Date.now();
const {browser, page} = await openStand();
try {
  for (const check of [checkRegistration, checkRoute, checkBrowse, checkEntitySource,
    checkDesignRowCap, checkViewerRowCap, checkDroppedFile, checkDroppedQuery]) {
    try {
      await check(page);
    } catch (e) {
      ok(check.name, false, `check threw: ${String(e).slice(0, 300)}`);
      await shot(page, `stand-${check.name}-failure`);
    }
  }
} finally {
  await browser.close();
}

note('stand/target', STAND_URL);
const code = report('designer.stand');
console.log(`${((Date.now() - started) / 1000).toFixed(1)}s · ${results.length} checks · artifacts in ${ARTIFACTS}`);
process.exit(code);
