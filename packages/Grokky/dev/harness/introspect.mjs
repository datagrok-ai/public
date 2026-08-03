#!/usr/bin/env node
// One-off: evaluates an expression in a logged-in Datagrok page. Used to settle "is this assert
// wrong or is the model wrong?" against the real API instead of by reading source.
//
//   node introspect.mjs "Array.from(DG.Viewer.fromType('Box plot', t).props.getProperties().map(p=>p.name))"
//
// `t` is Demographics; pass --table chem/SPGI or --table beer to use another demo file.

import * as fs from 'node:fs';
import * as path from 'node:path';
import {chromium} from 'playwright';
import yaml from 'yaml';

const expr = process.argv.filter((a) => !a.startsWith('--'))[2];
const tableArg = process.argv.indexOf('--table');
const table = tableArg >= 0 ? process.argv[tableArg + 1] : 'demog';
if (!expr) {
  console.error('usage: node introspect.mjs "<expression>" [--table demog|beer|chem/SPGI]');
  process.exit(1);
}

const cfg = yaml.parse(fs.readFileSync(
  path.join(process.env.USERPROFILE ?? process.env.HOME, '.grok', 'config.yaml'), 'utf8'));
const server = cfg.servers[cfg.default];
const apiUrl = server.url.replace(/\/$/, '');
const {token} = await (await fetch(`${apiUrl}/users/login/dev/${server.key}`, {method: 'POST'})).json();
const settings = await (await fetch(`${apiUrl}/admin/plugins/admin/settings`, {headers: {Authorization: token}})).json();
const webUrl = settings.settings.webRoot.replace(/\/$/, '');

const browser = await chromium.launch();
try {
  const ctx = await browser.newContext();
  const page = await ctx.newPage();
  await page.goto(webUrl + '/oauth/', {waitUntil: 'domcontentloaded', timeout: 120000});
  await ctx.addCookies([{name: 'auth', value: token, domain: new URL(webUrl).hostname, path: '/'}]);
  await page.evaluate((tk) => window.localStorage.setItem('auth', tk), token);
  await page.goto(webUrl, {waitUntil: 'domcontentloaded', timeout: 120000});
  await page.waitForFunction(() => document.querySelector('#grok-preloader, .grok-preloader') == null,
    null, {timeout: 120000}).catch(() => {});
  // grok.shell existing is not enough — the Dart interop table fills in progressively, and a
  // dapi call made too early throws "api.grok_Dapi_... is not a function". Wait on the actual
  // binding this script needs rather than on a proxy for it.
  await page.waitForFunction(() => typeof window.grok?.dapi?.files?.readAsText === 'function' &&
    typeof window.DG?.Viewer?.fromType === 'function', null, {timeout: 120000});

  const out = await page.evaluate(async ({expr, table}) => {
    const grok = window.grok; const DG = window.DG;
    const df = table === 'randomWalk' ? grok.data.demo.randomWalk() :
      await grok.dapi.files.readCsv(`System:DemoFiles/${table}.csv`);
    const view = grok.shell.addTableView(df);
    const t = view.dataFrame;
    await new Promise((r) => setTimeout(r, 300));
    try {
      // eslint-disable-next-line no-new-func
      const fn = new Function('grok', 'DG', 'view', 't', `return (async () => (${expr}))();`);
      const v = await fn(grok, DG, view, t);
      return {ok: true, value: JSON.stringify(v, null, 2)?.slice(0, 4000) ?? String(v)};
    } catch (e) {
      return {ok: false, value: String(e?.message ?? e)};
    } finally {
      try {
        view.close();
      } catch { /* view already gone */ }
    }
  }, {expr, table});

  console.log(out.ok ? out.value : `ERROR: ${out.value}`);
  process.exit(out.ok ? 0 : 1);
} finally {
  await browser.close();
}
