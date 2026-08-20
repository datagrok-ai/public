/* Lane 2 — what only a stand can answer: that the package as published (not as staged) registers
   the app, that the platform's own routes reach it, and that it is filed where Browse expects it.
   Everything about the designer's behaviour is Lane 1's job — do not grow this file into a copy of
   designer.spec.mjs. Run with `npm run e2e:stand` on a stand nobody else is using; publish first
   (`npm run debug-u2demo-<host>` in packages/U2Demo). */
import {ARTIFACTS, note, ok, report, results, shot} from './local.mjs';
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

const started = Date.now();
const {browser, page} = await openStand();
try {
  for (const check of [checkRegistration, checkRoute, checkBrowse]) {
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
