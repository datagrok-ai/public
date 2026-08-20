/* `npm run e2e:local` — preconditions, one browser, every check, one results.json.
   Preconditions are the whole difference between this lane and a stand: pub serve must be up
   (it is the only server involved) and the package must be staged into web/local/pkg/. */
import {execFileSync} from 'child_process';
import {existsSync, readdirSync} from 'fs';
import {join, resolve} from 'path';
import {ARTIFACTS, E2E_DIR, LOCAL_URL, note, ok, openApp, openLocal, report, results, shot} from './local.mjs';
import {APP, checks} from './designer.spec.mjs';

const REPO = resolve(E2E_DIR, '../../../..');
const PKG_DIR = join(REPO, 'core/client/xamgle/web/local/pkg');
const GROK_CORE = join(REPO, 'core/tools/grok-core/bin/grok-core.js');

async function pubServeIsUp() {
  try {
    const response = await fetch(`${LOCAL_URL}/login.html`, {signal: AbortSignal.timeout(15000)});
    return response.ok;
  } catch (e) {
    return false;
  }
}

/** Stages the package the first time; a rebuild is re-staged by `grok-core local watch`. */
function ensureStaged(name) {
  if (!existsSync(PKG_DIR) || !existsSync(GROK_CORE))
    return note('precondition/staging', `no local-mode tree at ${PKG_DIR} — assuming the shell is served elsewhere`);
  if (readdirSync(PKG_DIR).some((d) => d.toLowerCase() === name.toLowerCase()))
    return;
  console.log(`staging ${name} into local mode ...`);
  execFileSync(process.execPath, [GROK_CORE, 'local', 'stage', name], {stdio: 'inherit'});
}

const started = Date.now();
if (!await pubServeIsUp()) {
  console.error(`u2 e2e: nothing answers at ${LOCAL_URL} — start it with \`grok-core local up\` ` +
    '(one pub serve at a time: never start a second).');
  process.exit(2);
}
ensureStaged(APP.package);

const {browser, page} = await openLocal();
const booted = Date.now();
try {
  const via = await openApp(page, APP.package, APP.func);
  if (via === 'module') {
    note('precondition/app-open', `${APP.package}:${APP.func} is missing from the local fixture — ` +
      'opened through the package module; refresh with `grok-core local fixture --packages ' + APP.package + '`');
  }
  await page.waitForSelector('.u2-designer', {timeout: 60000});
  for (const check of checks) {
    try {
      await check.run(page);
    } catch (e) {
      ok(check.id, false, `check threw: ${String(e).slice(0, 300)}`);
      await shot(page, `${check.id.split(' ')[0]}-failure`);
    }
  }
} finally {
  await browser.close();
}

const code = report('designer.local');
console.log(`boot ${((booted - started) / 1000).toFixed(1)}s · run ${((Date.now() - booted) / 1000).toFixed(1)}s · ` +
  `${results.length} checks · artifacts in ${ARTIFACTS}`);
process.exit(code);
