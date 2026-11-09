/* The engine names three ownership mistakes that otherwise pass in silence: a control's field that takes the name of
   one of the control's own members (`name`, `root`, `effect`, …; strict mode), a child control held in a field but
   built outside the control's scope (a bare `input = new TextInput()` initializer; u2 debug on), and a control built
   under a scope that is already disposed. */

import {after, before, test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {debug, Usage, U2UsageError} from '../src/core/scope.js';
import {Scope, Control, Dialog, TextInput, AsyncValue, component, div, span, signal} from '../src/index.js';
import {Stream} from './platform-doubles.mjs';

before(() => debug.enable({stacks: false, orphanSeconds: 0}));
after(() => debug.disable());

async function consoleOf(body) {
  const said = {warn: [], error: []};
  const {warn, error} = console;
  console.warn = (...a) => said.warn.push(a.join(' '));
  console.error = (...a) => said.error.push(a.join(' '));
  try {
    await body();
    await flush();
  } finally {
    console.warn = warn;
    console.error = error;
    resetDom();
  }
  return said;
}

/** Runs `body` with the mode `mode`, and answers what it answers. */
function inMode(mode, body) {
  const was = Usage.mode;
  Usage.mode = mode;
  try {
    return body();
  } finally {
    Usage.mode = was;
  }
}

/** The misuse `body` makes, thrown: the mode is `throw` while it runs. */
function thrownBy(body) {
  try {
    inMode('throw', body);
  } catch (e) {
    return e;
  }
  assert.fail('no usage error was thrown');
}

/** What `body` reports in `report` mode, from a clean slate. */
function reportedBy(body) {
  const was = Usage.mode;
  Usage.mode = 'report';
  Usage.clear();
  try {
    body();
    return Usage.reports();
  } finally {
    Usage.mode = was;
  }
}

class NamedCard extends Control {
  name = 'Aspirin';
  effect = signal(0);
}

class RootCard extends Control {
  root = div([]);
}

/** What `body` reports in `report` mode once the microtask after it has run, from a clean slate. */
async function reportedAfter(body) {
  const was = Usage.mode;
  Usage.mode = 'report';
  Usage.clear();
  try {
    await body();
    await flush();
    return Usage.reports();
  } finally {
    Usage.mode = was;
    resetDom();
  }
}

test('reserved-key: a field named like one of the control\'s own members is reported, naming both', async () => {
  const live = Scope.liveCount;
  const reports = await reportedAfter(() => {
    new NamedCard().dispose();
    new RootCard().dispose();
  });
  const messages = reports.filter((r) => r.code === 'reserved-key').map((r) => r.message);
  assert.equal(messages.length, 3, messages.join('\n'));
  assert.match(messages[0], /NamedCard\.name overrides Control's own `name`/);
  assert.match(messages[1], /NamedCard\.effect overrides Control's own `effect\(\)`/);
  assert.match(messages[2], /RootCard\.root overrides Control's own `root`/);
  assert.equal(Scope.liveCount, live);
});

class FieldEditor extends Control {
  input = new TextInput({label: 'Name'});
  owned = this.build(() => new TextInput({label: 'Owned'}));

  constructor(given) {
    super();
    this.given = given;
    this.root.append(this.input.root, this.owned.root);
  }
}

test('a child control held in a field but built outside the control\'s scope warns, naming the field', async () => {
  const live = Scope.liveCount;
  let page;
  const said = await consoleOf(() => {
    const given = component(() => [new TextInput({label: 'Given'})]);
    page = component(() => div([new FieldEditor(given)]));
    given.dispose();
  });
  page.dispose();
  assert.equal(said.warn.length, 1, said.warn.join('\n'));
  assert.match(said.warn[0], /FieldEditor\.input holds a control built outside FieldEditor's scope/);
  assert.match(said.warn[0], /this\.build\(\(\) => /);
  assert.equal(Scope.liveCount, live);
});

class Holder extends Control {
  attach(x) {
    this.attached = x;
  }
}

test('a sibling handed to a control after its constructor is not its child: no warning', async () => {
  let page;
  const said = await consoleOf(() => {
    page = component(() => {
      const holder = new Holder();
      const sibling = new TextInput({label: 'Sibling'});
      holder.attach(sibling);
      return div([holder, sibling]);
    });
  });
  page.dispose();
  assert.deepEqual(said.warn, []);
});

test('with debug off nothing is checked', async () => {
  debug.disable();
  try {
    let said;
    await reportedAfter(async () => said = await consoleOf(() => {
      new NamedCard().dispose();
      component(() => div([new FieldEditor(undefined)])).dispose();
    }));
    assert.deepEqual(said, {warn: [], error: []});
  } finally {
    debug.enable({stacks: false, orphanSeconds: 0});
  }
});

test('usage: a misuse is reported once per call site, with its count, and logged nowhere on the console', async () => {
  const gone = Scope.root('session');
  gone.dispose();
  let reports;
  const said = await consoleOf(() => {
    reports = reportedBy(() => {
      for (let i = 0; i < 3; i++)
        gone.effect(() => {});
      gone.listen(document, 'click', () => {});
    });
  });
  assert.deepEqual(reports.map((r) => [r.code, r.count]), [['owner-disposed', 3], ['owner-disposed', 1]]);
  assert.match(reports[0].site, /dev-checks\.test\.js:\d+/);
  assert.notEqual(reports[0].site, reports[1].site);
  assert.deepEqual(said, {warn: [], error: []}, 'Node reports on stderr, past the console');
});

test('usage: off says nothing, throw throws a U2UsageError carrying the code', () => {
  const gone = Scope.root('session');
  gone.dispose();
  const was = Usage.mode;
  Usage.mode = 'off';
  Usage.clear();
  try {
    gone.setTimeout(() => {});
    assert.deepEqual(Usage.reports(), []);
  } finally {
    Usage.mode = was;
  }
  const e = thrownBy(() => gone.setInterval(() => {}, 10));
  assert.ok(e instanceof U2UsageError);
  assert.equal(e.code, 'owner-disposed');
  assert.match(e.message, /^u2 usage: owner-disposed: setInterval\(\) on a scope already disposed \(session\)/);
});

test('owner-missing: a control built with no owner reports once it owns something, naming the fixes', () => {
  const e = thrownBy(() => new TextInput({label: 'Loose'}));
  assert.equal(e.code, 'owner-missing');
  assert.match(e.message, /TextInput was built with no owner/);
  assert.match(e.message, /owner\.build\(\(\) => new TextInput\(…\)\)/);
  assert.match(e.message, /owner\.task\(async \(t\) => … t\.build\(…\)\)/);
  assert.match(e.message, /Scope\.root\('<what>'\)/);
});

test('owner-missing: none for a static control, a root, a built control, or one under Scope.root', () => {
  const root = Scope.root('page');
  assert.doesNotThrow(() => inMode('throw', () => {
    new Control().dispose();
    new Dialog('Root').dispose();
    component(() => [new TextInput({label: 'Owned'})]).dispose();
    root.build(() => new TextInput({label: 'Rooted'}));
  }));
  root.dispose();
});

test('owner-missing: an owned thing given no owner while nothing builds throws; a binding throws in any mode', () => {
  const e = thrownBy(() => new AsyncValue(async () => 1, {auto: false}));
  assert.equal(e.code, 'owner-missing');
  assert.match(e.message, /pass `owner`, or build it under `owner\.build\(…\)`/);
  const binding = (x) => x instanceof U2UsageError && x.code === 'owner-missing' &&
    /a signal binding needs an owner/.test(x.message) && /component\(\(self\) => …\)/.test(x.message);
  assert.throws(() => inMode('throw', () => span(signal('x'))), binding);
  assert.throws(() => inMode('report', () => span(signal('x'))), binding);
});

test('owner-disposed: a control built under a disposed owner, and every verb on one, names the fix', () => {
  const gone = Scope.root('session');
  gone.dispose();
  const built = thrownBy(() => gone.build(() => new TextInput({label: 'Late'})));
  assert.equal(built.code, 'owner-disposed');
  assert.match(built.message, /TextInput was built under a scope already disposed \(session\)/);
  assert.match(built.message, /Scope\.root\('<what>'\)/);
  const verbs = {
    effect: () => gone.effect(() => {}),
    listen: () => gone.listen(document, 'click', () => {}),
    setTimeout: () => gone.setTimeout(() => {}),
    requestFrame: () => gone.requestFrame(() => {}),
    subscribe: () => gone.subscribe(new Stream(), () => {}),
    task: () => gone.task(() => {}),
  };
  for (const [what, call] of Object.entries(verbs)) {
    const e = thrownBy(call);
    assert.equal(e.code, 'owner-disposed', what);
    assert.ok(e.message.includes(`${what}() on a scope already disposed (session)`), e.message);
    assert.match(e.message, /owner\.task\(async \(t\) => … t\.build\(…\)\)/);
  }
});

test('a subclass of a library class is not reported; a reserved field is, in strict mode with debug off', async () => {
  class MyInput extends TextInput {}
  class CampaignInput extends TextInput {
    name = 'campaign';
  }
  debug.disable();
  const owner = Scope.root('page');
  try {
    const reports = await reportedAfter(() => owner.build(() => {
      new MyInput({label: 'Mine'});
      new CampaignInput({label: 'Campaign'});
    }));
    assert.deepEqual(reports.map((r) => r.code), ['reserved-key']);
    assert.match(reports[0].message, /^CampaignInput\.name overrides Control's own `name`/);
  } finally {
    owner.dispose();
    debug.enable({stacks: false, orphanSeconds: 0});
  }
});

test('reserved-key: a component() handle key that is a Control member or an owner verb throws; in production it ' +
  'is not copied', () => {
  for (const [key, handle] of [['name', {name: 'mine'}], ['dispose', {dispose: () => {}}], ['build', {build: 1}]]) {
    const e = thrownBy(() => component(() => ({root: div([]), ...handle})));
    assert.equal(e.code, 'reserved-key');
    assert.ok(e.message.includes(`component()'s handle has \`${key}\`, a member every control has`), e.message);
  }
  const built = inMode('report', () => component(() => ({root: div([]), dispose: () => {}, build: 1, count: 1})));
  assert.equal(built.count, 1);
  assert.equal(built.dispose, Control.prototype.dispose, 'the control keeps its own member');
  assert.equal(typeof built.build, 'function', 'and its owner verb');
  built.dispose();
});

test('source-missing: subscribing to an event the host does not have throws; in production it subscribes nothing', () => {
  const owner = Scope.root('page');
  try {
    const e = thrownBy(() => owner.subscribe(undefined, () => {}));
    assert.equal(e.code, 'source-missing');
    assert.match(e.message, /^u2 usage: source-missing: owner\.subscribe\(source, handler\) needs an observable; got undefined/);
    assert.match(thrownBy(() => owner.subscribe({}, () => {})).message, /got an object with no subscribe method/);
    assert.equal(typeof inMode('report', () => owner.subscribe({}, () => {})), 'function');
  } finally {
    owner.dispose();
  }
});
