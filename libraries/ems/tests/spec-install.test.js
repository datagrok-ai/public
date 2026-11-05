import {emsBackends} from '../src/core/domain-backend.js';
import {test, testOwner} from '@datagrok-libraries/u2/tests/owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {readdirSync, readFileSync, existsSync} from 'node:fs';
import {fileURLToPath} from 'node:url';
import {join, sep} from 'node:path';
import {flush, resetDom} from '@datagrok-libraries/u2/tests/dom-shim.js';
import {backends, Scope, Registry, SpecContext, renderSpec, registerAll, notify} from '@datagrok-libraries/u2';
import {MemoryDomainBackend} from '../src/core/memory-domain.js';
import {SharedSession} from '../src/core/session.js';
import {emsSourcesPlugin} from '../src/core/registrations.js';
import {Property, WidgetDescriptor, platform} from '@datagrok-libraries/u2/tests/platform-doubles.mjs';
import {backend} from './domain-fixtures.mjs';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {registerPlatformComponents} = await import('@datagrok-libraries/u2/dg');
const {installEms, emsPlugin, registerEmsComponents} = await import('../src/dg/registrations.js');

const PACKAGES = fileURLToPath(new URL('../../../packages/', import.meta.url));
const RECIPES = fileURLToPath(new URL('../docs/recipes/', import.meta.url));

function scoped(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    const saved = {...backends};
    const savedEms = {...emsBackends};
    WidgetDescriptor.registry = [new WidgetDescriptor('Grid', [new Property('allowEdit', 'bool', {defaultValue: false})])];
    try {
      await body();
    } finally {
      for (const key of Object.keys(backends))
        delete backends[key];
      Object.assign(backends, saved);
      for (const key of Object.keys(emsBackends))
        delete emsBackends[key];
      Object.assign(emsBackends, savedEms);
      WidgetDescriptor.registry = [];
      platform.reset();
      notify.closeAll();
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

async function settle() {
  for (let i = 0; i < 5; i++)
    await flush();
}

const isPlaceholder = (built) => built.classList?.contains('u2-spec-error') === true;
const domainTags = (reg) => reg.metas().map((m) => m.tag).filter((t) => t.startsWith('u2-domain-') || t.startsWith('u2-manifest-'));

scoped('installEms: the source, every control and the session in one install — once, and undone by its uninstall', async () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  assert.deepEqual(domainTags(reg), [], 'u2 registers no domain tag');
  const uninstall = installEms(reg);
  const tags = domainTags(reg);
  assert.deepEqual(['u2-domain-source', 'u2-domain-form', 'u2-domain-list', 'u2-domain-pick', 'u2-domain-grid',
    'u2-domain-data-table', 'u2-domain-tree', 'u2-domain-history', 'u2-domain-children', 'u2-domain-search',
    'u2-domain-filters', 'u2-manifest-editor'].filter((t) => !tags.includes(t)), []);
  assert.deepEqual(reg.plugins(), [emsPlugin]);
  assert.equal(installEms(reg), uninstall, 'a second install is the first');
  registerEmsComponents(reg);
  assert.deepEqual(reg.plugins(), [emsPlugin], 'the deprecated alias installs the same plugin');
  uninstall();
  assert.deepEqual(domainTags(reg), []);
  assert.deepEqual(reg.plugins(), []);
});

scoped('P94: two sources of one spec share one session after an unrelated use() — and without ems none', async () => {
  emsBackends.domain = backend();
  const reg = new Registry();
  registerAll(reg);
  installEms(reg);
  let wrapped = 0;
  reg.use({name: 'unrelated', register: () => {}, specAmbient: (instance, build) => {
    wrapped++;
    return build();
  }});
  const spec = {$schema: 'dg-ui/1', root: {tag: 'div'}, components: [
    {tag: 'u2-domain-source', name: 'projects', props: {table: 'grit.project'}},
    {tag: 'u2-domain-source', name: 'issues', props: {table: 'grit.issue'}},
  ]};
  const instance = renderSpec(spec, new SpecContext(), reg);
  await settle();
  const session = SharedSession.of(instance);
  assert.ok(session instanceof SharedSession);
  assert.equal(instance.node('projects').session, session);
  assert.equal(instance.node('issues').session, session, 'one Save, one transaction');
  assert.equal(wrapped, 2, 'the other plugin wraps every source too');
  assert.equal(SharedSession.ambient, undefined);
  instance.dispose();

  const headless = new Registry();
  registerAll(headless);
  headless.use(emsSourcesPlugin);
  const alone = testOwner().build(() => renderSpec(spec, new SpecContext(), headless));
  await settle();
  assert.equal(alone.node('issues').session, SharedSession.of(alone), 'the headless install brings the session');
  alone.dispose();
});

/** Every spec the repo ships: the packages' `*.spec.json`, and the ```json blocks of the recipes. */
function shippedSpecs() {
  const specs = [];
  const walk = (dir) => {
    for (const entry of readdirSync(dir, {withFileTypes: true})) {
      if (entry.isDirectory() && entry.name !== 'node_modules')
        walk(join(dir, entry.name));
      else if (entry.name.endsWith('.spec.json'))
        specs.push({from: join(dir, entry.name), spec: JSON.parse(readFileSync(join(dir, entry.name), 'utf8'))});
    }
  };
  for (const pkg of readdirSync(PACKAGES)) {
    if (existsSync(join(PACKAGES, pkg, 'src')))
      walk(join(PACKAGES, pkg, 'src'));
  }
  for (const file of readdirSync(RECIPES).filter((f) => f.endsWith('.md'))) {
    const text = readFileSync(join(RECIPES, file), 'utf8');
    for (const [, block] of text.matchAll(/```json\r?\n([\s\S]*?)```/g)) {
      const spec = JSON.parse(block);
      if (spec.$schema === 'dg-ui/1')
        specs.push({from: `docs/recipes/${file}`, spec});
    }
  }
  return specs;
}

/** The schema.json a package ships for the tables a spec's sources name. */
function schemaFor(spec) {
  const names = new Set((spec.components ?? []).filter((n) => n.tag === 'u2-domain-source')
    .map((n) => String(n.props?.table ?? '').split('.')[0]));
  assert.ok(names.size <= 1, `one schema per spec: ${[...names].join(', ')}`);
  const [name] = names;
  if (name === undefined)
    return undefined;
  for (const pkg of readdirSync(PACKAGES)) {
    const path = join(PACKAGES, pkg, 'databases', name, 'schema.json');
    if (existsSync(path))
      return JSON.parse(readFileSync(path, 'utf8'));
  }
  assert.fail(`no package ships the "${name}" schema`);
}

const SPECS = shippedSpecs();

test('the shipped specs are found: Stockroom\'s app and the spec-app recipe', () => {
  const from = SPECS.map((s) => s.from);
  assert.ok(from.some((f) => f.split(sep).join('/').endsWith('Stockroom/src/app.spec.json')), from.join(', '));
  assert.ok(from.includes('docs/recipes/spec-app.md'), from.join(', '));
});

for (const {from, spec} of SPECS) {
  scoped(`P6: ${from.replace(PACKAGES, 'packages/')} renders on the memory backend with no placeholder`, async () => {
    const schema = schemaFor(spec);
    if (schema !== undefined)
      emsBackends.domain = new MemoryDomainBackend(schema);
    const reg = new Registry();
    registerPlatformComponents(reg);
    installEms(reg);
    const warnings = [];
    const original = console.warn;
    console.warn = (message) => warnings.push(String(message));
    let instance;
    try {
      instance = renderSpec(spec, new SpecContext(), reg);
      document.body.append(instance.root);
      await settle();
    } finally {
      console.warn = original;
    }
    const broken = [...instance.nodes()].filter(([, built]) => isPlaceholder(built))
      .map(([node, built]) => `${node.name ?? node.tag}: ${built.textContent}`);
    assert.deepEqual(broken, []);
    assert.deepEqual(warnings.filter((w) => w.startsWith('u2 spec:')), []);
    if ((spec.components ?? []).some((n) => n.tag === 'u2-domain-source'))
      assert.ok(SharedSession.of(instance) instanceof SharedSession, 'the sources joined the spec session');
    instance.dispose();
  });
}
