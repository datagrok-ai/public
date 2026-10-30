// u2 registers no domain tag; registerEmsComponents adds them

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {Registry, registerAll} from '@datagrok-libraries/u2';
import {registerDomainSources} from '../src/core/registrations.js';

register('@datagrok-libraries/u2/tests/dg-stub.mjs', import.meta.url);
const {registerPlatformComponents} = await import('@datagrok-libraries/u2/dg');
const {registerEmsComponents} = await import('../src/dg/registrations.js');

const CONTROLS = ['u2-domain-form', 'u2-domain-list', 'u2-domain-pick', 'u2-domain-grid',
  'u2-domain-data-table', 'u2-domain-tree', 'u2-domain-history', 'u2-domain-children',
  'u2-domain-search', 'u2-domain-filters'];

const tagsOf = (reg) => reg.metas().map((m) => m.tag);

test('neither u2 entry point registers a domain tag', () => {
  const core = new Registry();
  registerAll(core);
  assert.deepEqual(tagsOf(core).filter((t) => t.startsWith('u2-domain-')), []);
  const platform = new Registry();
  registerPlatformComponents(platform);
  assert.deepEqual(tagsOf(platform).filter((t) => t.startsWith('u2-domain-')), []);
});

test('registerDomainSources adds the source alone — what a headless consumer needs', () => {
  const reg = new Registry();
  registerDomainSources(reg);
  assert.deepEqual(tagsOf(reg), ['u2-domain-source']);
  const meta = reg.get('u2-domain-source');
  assert.equal(meta.category, 'Data');
  assert.equal(meta.visual, false);
});

test('registerEmsComponents adds the source and every control, and is idempotent', () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  registerEmsComponents(reg);
  const tags = tagsOf(reg);
  assert.deepEqual(['u2-domain-source', ...CONTROLS].filter((t) => !tags.includes(t)), []);
  registerEmsComponents(reg);
  assert.deepEqual(tagsOf(reg), tags, 'a second call leaves the registry as it was');
});
