// Evaluates the built js-api bundle in a Node sandbox with stubbed browser globals and returns
// its `DG` export. Shared by the smoke test and the unit tests under scripts/unit.
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const DEFAULT_BUNDLE = path.resolve(__dirname, '../../../core/client/xamgle/web/js/api/js-api.js');

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

/** Loads the bundle at [bundlePath] (default: the checked-in xamgle copy). Throws an Error whose
 * `frames` property holds the first call frames when the bundle fails to initialize. */
function loadBundle(bundlePath = DEFAULT_BUNDLE) {
  if (!fs.existsSync(bundlePath))
    throw new Error('bundle not found - run the webpack build first: ' + bundlePath);

  // Evaluating against stubs makes deferred work (timers, promises) fail in ways a real
  // browser never would. The synchronous module-init phase is the only part that matters.
  process.on('uncaughtException', () => {});
  process.on('unhandledRejection', () => {});

  const sandbox = {console, setTimeout, clearTimeout, setInterval, clearInterval};
  sandbox.globalThis = sandbox;
  sandbox.window = sandbox;
  sandbox.self = sandbox;
  // Deferred init code registers window listeners; give it real no-ops so it cannot throw
  // after a test file has finished.
  sandbox.addEventListener = () => {};
  sandbox.removeEventListener = () => {};
  sandbox.dispatchEvent = () => true;
  sandbox.requestAnimationFrame = (f) => setTimeout(f, 0);
  for (const name of ['document', 'navigator', 'location', 'rxjs', 'OCL', 'fetch', 'localStorage'])
    sandbox[name] = stub(name);
  vm.createContext(sandbox);

  try {
    new vm.Script(fs.readFileSync(bundlePath, 'utf8'), {filename: 'js-api.js'})
      .runInContext(sandbox, {timeout: 120000});
  } catch (e) {
    // The bundle is one minified line, so V8 splices the whole source into the stack as
    // context. Keep only the call frames — their offsets locate the throwing class.
    const err = new Error('the bundle threw while initializing: ' + e.name + ': ' + e.message);
    err.frames = (e.stack || '').split('\n').filter((line) => /^\s+at /.test(line)).slice(0, 5).join('\n');
    throw err;
  }
  return {DG: sandbox.DG, sandbox, bundlePath};
}

module.exports = {loadBundle, DEFAULT_BUNDLE};
