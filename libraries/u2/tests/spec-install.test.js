/* Plugins (WO-C14, P92, P94): a companion library installs into a registry with `use` — its tags
   and its spec ambient together, once per name, undone by the uninstall it answers — and a spec
   naming a tag of a library the host has not installed says which one to install. The "ems" here
   is a stand-in: u2 never imports the real one. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Component} from '../src/core/component.js';
import {Registry} from '../src/spec/registry.js';
import {SpecContext, renderSpec} from '../src/spec/spec.js';
import {registerAll} from '../src/spec/registrations.js';
import * as u2 from '../src/index.js';

function install(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

function captureWarnings(body) {
  const warnings = [];
  const original = console.warn;
  console.warn = (message) => warnings.push(String(message));
  try {
    return {result: body(), warnings};
  } finally {
    console.warn = original;
  }
}

/** A tray tag whose component records what it was built under. */
function trayTag(tag, built) {
  return {tag, visual: false, description: tag, props: [], example: {tag, name: 'x'},
    createComponent: () => {
      const component = new Component();
      component.ambient = [...stack];
      component.session = currentSession;
      built.push(component);
      return component;
    }};
}

/** The ambients that are open around a build, outermost first. */
const stack = [];
let currentSession;

function plugin(name, tags, extra = {}) {
  return {name, register: (reg) => {
    for (const meta of tags) {
      if (reg.get(meta.tag) === undefined)
        reg.register(meta);
    }
  }, specAmbient: (instance, build) => {
    stack.push(name);
    try {
      return build();
    } finally {
      stack.pop();
    }
  }, ...extra};
}

/** What ems does: one session per spec instance, ambient while a tray component builds. */
function fakeEms(built) {
  const sessions = new WeakMap();
  return {
    sessions,
    plugin: {
      name: '@datagrok-libraries/ems',
      tagPrefixes: ['u2-domain-', 'u2-manifest-'],
      register: (reg) => reg.register(trayTag('u2-domain-source', built)),
      specAmbient: (instance, build) => {
        let session = sessions.get(instance);
        if (session === undefined)
          sessions.set(instance, session = {instance});
        const outer = currentSession;
        currentSession = session;
        try {
          return build();
        } finally {
          currentSession = outer;
        }
      },
    },
  };
}

const spec = (components, root = {tag: 'div'}) => ({$schema: 'dg-ui/1', components, root});

install('two plugins\' ambients nest in install order, the first installed outermost', () => {
  const built = [];
  const reg = new Registry();
  reg.use(plugin('outer', [trayTag('u2-outer-source', built)]));
  reg.use(plugin('inner', [trayTag('u2-inner-source', built)]));
  assert.deepEqual(reg.plugins().map((p) => p.name), ['outer', 'inner']);
  const instance = renderSpec(spec([{tag: 'u2-outer-source', name: 'a'}, {tag: 'u2-inner-source', name: 'b'}]),
    new SpecContext(), reg);
  assert.deepEqual(built.map((c) => c.ambient), [['outer', 'inner'], ['outer', 'inner']],
    'every tray build runs under every plugin\'s ambient, whichever plugin owns the tag');
  assert.deepEqual(stack, [], 'nothing stays open past the build');
  instance.dispose();
});

install('uninstall removes the tags the plugin registered and its ambient — and nothing else', () => {
  const built = [];
  const reg = new Registry();
  registerAll(reg);
  const before = reg.metas().length;
  const shared = trayTag('u2-shared-source', built);
  reg.register(shared);
  const uninstall = reg.use(plugin('lib', [trayTag('u2-lib-source', built), shared]));
  assert.notEqual(reg.get('u2-lib-source'), undefined);
  assert.equal(reg.use(reg.plugins()[0]), uninstall, 'installed already: the same uninstall');
  uninstall();
  assert.equal(reg.get('u2-lib-source'), undefined, 'its own tag is gone');
  assert.notEqual(reg.get('u2-shared-source'), undefined, 'a tag registered before the install stays');
  assert.equal(reg.metas().length, before + 1);
  assert.deepEqual(reg.plugins(), []);
  uninstall();
  const instance = renderSpec(spec([{tag: 'u2-shared-source', name: 'a'}]), new SpecContext(), reg);
  assert.deepEqual(built.map((c) => c.ambient), [[]], 'the ambient is unhooked');
  instance.dispose();
});

install('a second copy of a plugin joins the first: one install, one warning, the first one\'s uninstall', () => {
  const reg = new Registry();
  const first = reg.use(plugin('lib', [trayTag('u2-lib-source', [])]));
  const {result: second, warnings} = captureWarnings(() => reg.use(plugin('lib', [trayTag('u2-lib-other', [])])));
  assert.equal(second, first);
  assert.equal(reg.get('u2-lib-other'), undefined, 'the second copy registers nothing');
  assert.equal(reg.plugins().length, 1);
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /"lib" is installed already/);
  first();
});

install('a domain tag without ems is a placeholder that names ems — in the form and on the tray', async () => {
  const reg = new Registry();
  registerAll(reg);
  const {result: instance, warnings} = captureWarnings(() => renderSpec(spec(
    [{tag: 'u2-domain-source', name: 'substances', props: {table: 'stockroom.substance'}}],
    {tag: 'u2-panel', children: [{tag: 'u2-domain-list', name: 'list'}, {tag: 'u2-manifest-editor', name: 'editor'}]}),
  new SpecContext(), reg));
  await flush();
  const placeholders = [...instance.root.querySelectorAll('.u2-spec-error')].map((el) => el.textContent);
  assert.equal(placeholders.length, 2);
  for (const text of placeholders)
    assert.match(text, /@datagrok-libraries\/ems — call installEms\(\)/);
  const tray = instance.node('substances');
  assert.match(tray.textContent, /u2-domain-source: is not a registered component; it comes with @datagrok-libraries\/ems/);
  assert.equal(warnings.some((w) => /u2-domain-source: .*installEms\(\)/.test(w)), true,
    'a tray placeholder never reaches the canvas, so the console says it');
  instance.dispose();
});

install('a tag the installed plugin does not have names the plugin, not the install', async () => {
  const reg = new Registry();
  registerAll(reg);
  reg.use(fakeEms([]).plugin);
  const {result: instance} = captureWarnings(() =>
    renderSpec(spec([], {tag: 'u2-domain-lst', name: 'typo'}), new SpecContext(), reg));
  await flush();
  assert.match(instance.root.querySelector('.u2-spec-error').textContent,
    /u2-domain-lst: .*@datagrok-libraries\/ems is installed and registers no such tag/);
  instance.dispose();
});

install('P94: two u2-domain-source nodes share one session after an unrelated use()', () => {
  const built = [];
  const ems = fakeEms(built);
  const reg = new Registry();
  registerAll(reg);
  reg.use(ems.plugin);
  reg.use(plugin('unrelated', []));
  const instance = renderSpec(spec([
    {tag: 'u2-domain-source', name: 'projects'},
    {tag: 'u2-domain-source', name: 'issues'},
  ]), new SpecContext(), reg);
  const session = ems.sessions.get(instance);
  assert.notEqual(session, undefined);
  assert.equal(instance.node('projects').session, session);
  assert.equal(instance.node('issues').session, session, 'one unit of work');
  assert.deepEqual(instance.node('issues').ambient, ['unrelated'], 'the other plugin wraps too');
  assert.equal(currentSession, undefined, 'nothing leaks past the build');

  const other = renderSpec(spec([{tag: 'u2-domain-source', name: 'more'}]), new SpecContext(), reg);
  assert.notEqual(other.node('more').session, session, 'another instance, another session');
  other.dispose();
  instance.dispose();
});

test('the single ambient slot is gone from the package root; the plugin type is the way in', () => {
  assert.equal('setSpecAmbient' in u2, false);
  assert.equal(typeof u2.Registry.prototype.use, 'function');
});
