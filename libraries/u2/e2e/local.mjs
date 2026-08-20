/* Lane 1 — the real Datagrok client with no server: `?mode=local` boots from static files and a
   constant token, so there is no login step and no stand to wait for. Apps are opened by applying
   their package function directly, never by navigating the UI. Also the small reporting surface both
   lanes share (see stand.mjs). */
import {mkdirSync, writeFileSync} from 'fs';
import {dirname, join} from 'path';
import {fileURLToPath} from 'url';
import {chromium} from 'playwright';

export const E2E_DIR = dirname(fileURLToPath(import.meta.url));
export const ARTIFACTS = join(E2E_DIR, '.artifacts');
export const LOCAL_URL = process.env.U2_LOCAL_URL ?? 'http://localhost:63343';
/** Local mode has no deployment, so the only real waits are compile-on-demand and the boot itself. */
export const BOOT_TIMEOUT = Number(process.env.U2_BOOT_TIMEOUT ?? 300000);

export const results = [];
export const consoleErrors = [];
export const pageErrors = [];

export function ok(id, pass, detail = '') {
  results.push({id, pass: !!pass, detail: String(detail).slice(0, 400)});
  console.log(`${pass ? 'PASS' : 'FAIL'}  ${id}  ${String(detail).slice(0, 400)}`);
  return !!pass;
}

export function note(id, detail = '') {
  results.push({id, pass: 'note', detail: String(detail).slice(0, 400)});
  console.log(`NOTE  ${id}  ${String(detail).slice(0, 400)}`);
}

/** Resource and rendering chatter a healthy client emits anyway — anything else is a real error. */
const NOISE = /Failed to load resource|favicon|ERR_CONNECTION|drawMolecule|Stack trace|^\s*packages\//;

export function capture(page) {
  page.on('console', (m) => {
    const text = m.text();
    if (m.type() === 'error' && !NOISE.test(text))
      consoleErrors.push(text.slice(0, 300));
    if (m.type() === 'warning' && /command|cmd:/i.test(text))
      consoleErrors.push('[warn] ' + text.slice(0, 300));
  });
  page.on('pageerror', (e) => {
    pageErrors.push(String(e).slice(0, 400));
    console.log('  [pageerror]', String(e).slice(0, 300));
  });
}

export async function launch({headless = process.env.U2_HEADED !== '1'} = {}) {
  const browser = await chromium.launch({channel: 'chrome', headless});
  const ctx = await browser.newContext({viewport: {width: 1600, height: 1000},
    permissions: ['clipboard-read', 'clipboard-write']});
  const page = await ctx.newPage();
  // the client is already up by the time checks run: a slow action means a broken one, not a busy one
  page.setDefaultTimeout(15000);
  capture(page);
  return {browser, page};
}

/** The shell is ready when the platform globals are up and the (constant-token) user is resolved. */
export async function openLocal(options = {}) {
  const {browser, page} = await launch(options);
  await page.goto(`${LOCAL_URL}/login.html?mode=local`, {waitUntil: 'load', timeout: BOOT_TIMEOUT});
  await page.waitForFunction(() => {
    try {
      return !!(window.DG && DG.Func && grok.shell.user);
    } catch (e) {
      return false;
    }
  }, null, {timeout: BOOT_TIMEOUT});
  return {browser, page};
}

/** Opens an app by applying its function — no Browse navigation, no routing, no publish. Function
 * metadata comes from the captured fixture, so a function added since the last
 * `grok-core local fixture` is unknown to the client: the package module is then loaded and the
 * export called directly. Answers which path was taken. */
export async function openApp(page, pkg, funcName) {
  return page.evaluate(async ({pkg, funcName}) => {
    const found = DG.Func.find({package: pkg, name: funcName});
    if (found.length) {
      grok.shell.addView(await found[0].apply());
      return 'func';
    }
    const sibling = DG.Func.find({package: pkg})[0];
    if (!sibling)
      throw new Error(`package ${pkg} is not in the local fixture — stage it and refresh the fixture`);
    await sibling.package.load();
    const module = sibling.package.getModule('package.js');
    if (typeof module?.[funcName] !== 'function')
      throw new Error(`${pkg}:${funcName} is neither registered nor exported by the staged bundle`);
    grok.shell.addView(await module[funcName]());
    return 'module';
  }, {pkg, funcName});
}

export function shot(page, name) {
  mkdirSync(ARTIFACTS, {recursive: true});
  return page.screenshot({path: join(ARTIFACTS, `${name}.png`)});
}

/** Evidence a check leaves behind whole — a `detail` string is truncated, a dump is not. */
export function artifact(name, text) {
  mkdirSync(ARTIFACTS, {recursive: true});
  writeFileSync(join(ARTIFACTS, name), text);
}

/** Writes the run to `.artifacts/<name>.json` and answers the process exit code. */
export function report(name = 'results') {
  mkdirSync(ARTIFACTS, {recursive: true});
  const failed = results.filter((r) => r.pass === false);
  writeFileSync(join(ARTIFACTS, `${name}.json`),
    JSON.stringify({results, consoleErrors, pageErrors}, null, 1));
  console.log(`\n${results.filter((r) => r.pass === true).length} passed, ${failed.length} failed, ` +
    `${consoleErrors.length} console errors, ${pageErrors.length} page errors`);
  for (const f of failed)
    console.log(`  FAIL ${f.id}: ${f.detail}`);
  for (const e of pageErrors)
    console.log(`  pageerror: ${e}`);
  return failed.length > 0 || pageErrors.length > 0 ? 1 : 0;
}
