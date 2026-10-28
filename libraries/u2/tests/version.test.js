import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {register} from 'node:module';
import {Scope, version} from 'datagrok-api/u2core';
import {VERSION} from '../src/index.js';

// `?engine=<v>` on a u2 module (and on everything it imports) serves it a u2core whose `version` is <v>, or missing
register('data:text/javascript,' + encodeURIComponent(`
const ENGINE = /[?&]engine=([^&]*)/;
export async function resolve(specifier, context, next) {
  const engine = context.parentURL?.startsWith('file:') ? ENGINE.exec(new URL(context.parentURL).search) : null;
  if (engine === null)
    return next(specifier, context);
  if (specifier === 'datagrok-api/u2core') {
    const real = await next(specifier, context);
    const v = engine[1] === 'none' ? 'undefined' : JSON.stringify(engine[1]);
    const source = 'export * from ' + JSON.stringify(real.url) + '; export const version = ' + v + ';';
    return {url: 'data:text/javascript,' + encodeURIComponent(source), format: 'module', shortCircuit: true};
  }
  if (specifier.startsWith('.')) {
    const url = new URL(specifier, context.parentURL);
    url.search = new URL(context.parentURL).search;
    return next(url.href, context);
  }
  return next(specifier, context);
}`));

const refusal = (v) => ({message: `u2 needs the u2 engine 1.2+; this platform has ${v}: update the platform`});

test('VERSION is the package.json version', () => {
  const pkg = JSON.parse(readFileSync(new URL('../package.json', import.meta.url), 'utf8'));
  assert.equal(VERSION, pkg.version);
});

test('the engine u2 is built on is 1.2 or newer', () => {
  const [major, minor] = version.split('.').map(Number);
  assert.ok(major > 1 || major === 1 && minor >= 2, version);
});

const header = (u2) => new u2.AppHeader({icon: 'flask', items: ['Browse the vault content.']});

test('on an engine without a version, u2 loads and building a component fails naming the fix', async () => {
  const u2 = await import('../src/index.js?engine=none');
  assert.throws(() => header(u2), refusal('an older one'));
  assert.throws(() => u2.component(() => u2.div()), refusal('an older one'));
  assert.throws(() => u2.Control.build(() => u2.div()), refusal('an older one'));
  assert.throws(() => new u2.TableSource(), refusal('an older one'));
});

test('on an engine older than 1.2 u2 loads and a control refuses naming its version; 1.2 and newer build', async () => {
  for (const v of ['1.1.0', '1.0.9', '0.9.0']) {
    const u2 = await import(`../src/index.js?engine=${v}`);
    assert.throws(() => header(u2), refusal(v));
  }
  for (const v of ['1.2.0', '1.12.0', '2.0.0']) {
    const u2 = await import(`../src/index.js?engine=${v}`);
    const owner = Scope.root('engine version');
    assert.ok(owner.build(() => header(u2)).root.classList.contains('u2-app-header'));
    owner.dispose();
  }
});
