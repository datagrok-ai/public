#!/usr/bin/env node
/**
 * Fails the build when the bundled js-api.js is dead on arrival.
 *
 * A bundle can download perfectly, be exactly the right size, and still throw the moment
 * a browser evaluates it. An import cycle that leaves a base class in the temporal dead
 * zone is enough: `class FormViewer extends Viewer` in a module emitted before viewer.ts
 * throws "Cannot access 'wi' before initialization" at init, window.grok is never
 * defined, and every user gets a stuck login screen. webpack does not warn about it,
 * tsc cannot see it, and HTTP health checks stay green because the asset itself is fine.
 *
 * Which side of a cycle webpack emits first is not stable — it shifts when unrelated
 * imports change — so this is not a one-off. The only reliable gate is to evaluate the
 * real bundle and require that it finishes and exports its core classes.
 *
 * Usage: node scripts/smoke-bundle.cjs [path/to/js-api.js]
 */
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const DEFAULT_BUNDLE = path.resolve(__dirname, '../../../core/client/xamgle/web/js/api/js-api.js');
const bundlePath = process.argv[2] ? path.resolve(process.argv[2]) : DEFAULT_BUNDLE;

// Evaluating against stubs makes deferred work (timers, promises) fail in ways a real
// browser never would. The synchronous module-init phase is the only part that matters
// and try/catch already covers it, so async noise is dropped.
process.on('uncaughtException', () => {});
process.on('unhandledRejection', () => {});

// Every property access yields another callable/constructible stub, so init code that
// pokes at DOM or rxjs globals keeps going instead of failing for an unrelated reason.
function stub(name) {
  const fn = function () { return stub(name + '()'); };
  return new Proxy(fn, {
    get(target, prop) {
      if (prop === 'then') return undefined;
      if (prop === 'prototype') return target.prototype;
      if (prop === Symbol.toPrimitive) return () => name;
      if (prop === Symbol.iterator) return function* () {};
      return stub(name + '.' + String(prop));
    },
    set: () => true,
    has: () => true,
    apply: () => stub(name + '()'),
    construct: () => stub('new ' + name),
  });
}

const REQUIRED = ['Viewer', 'Grid', 'FormViewer', 'Point', 'DataFrame', 'Column'];

function fail(message, detail) {
  console.error('\njs-api bundle smoke test FAILED');
  console.error('  bundle: ' + bundlePath);
  console.error('  ' + message);
  if (detail)
    console.error('\n' + detail);
  console.error('\nThis bundle is served to every browser as /js/api/js-api.js. In this state the');
  console.error('app cannot start: window.grok is never defined and users get a stuck login page.');
  console.error('A likely cause is a new import cycle whose base class is evaluated too early —');
  console.error('make the offending import `import type`, or reach the class through DG.* at');
  console.error('call time, so the runtime edge disappears.');
  process.exit(1);
}

if (!fs.existsSync(bundlePath))
  fail('bundle not found - run the webpack build first');

const sandbox = {console, setTimeout, clearTimeout, setInterval, clearInterval};
sandbox.globalThis = sandbox;
sandbox.window = sandbox;
sandbox.self = sandbox;
for (const name of ['document', 'navigator', 'location', 'rxjs', 'OCL', 'fetch', 'localStorage'])
  sandbox[name] = stub(name);
vm.createContext(sandbox);

try {
  new vm.Script(fs.readFileSync(bundlePath, 'utf8'), {filename: 'js-api.js'})
    .runInContext(sandbox, {timeout: 120000});
} catch (e) {
  // The bundle is one minified line, so V8 splices the whole source into the stack as
  // context. Keep only the call frames — their offsets locate the throwing class.
  const frames = (e.stack || '').split('\n')
    .filter((line) => /^\s+at /.test(line))
    .slice(0, 5)
    .join('\n');
  fail('the bundle threw while initializing: ' + e.name + ': ' + e.message, frames);
}

const DG = sandbox.DG;
if (!DG || typeof DG !== 'object')
  fail('the bundle evaluated but never defined the DG export');

const missing = REQUIRED.filter((name) => typeof DG[name] !== 'function');
if (missing.length)
  fail('DG is missing expected classes: ' + missing.join(', '));

if (typeof DG.U2?.Control !== 'function' || typeof DG.U2?.signal !== 'function')
  fail('DG.U2 is missing its u2core exports (Control, signal)');

for (const name of ['Widget', 'Viewer']) {
  if (typeof DG[name] !== 'function' || !(DG[name].prototype instanceof DG.U2.Control))
    fail('DG.' + name + ' does not extend DG.U2.Control');
}

console.log(`js-api bundle smoke test passed (${REQUIRED.length} core exports present)`);
