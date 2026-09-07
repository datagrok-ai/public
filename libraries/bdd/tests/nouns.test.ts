import assert from 'node:assert/strict';
import {beforeEach, test} from 'node:test';
import {contextFirst, describeNoun, NounError, parseNoun} from '../src/nouns.js';
import {context, element, kind, resetRegistry} from '../src/registry.js';

beforeEach(() => {
  resetRegistry();
  kind('button', {selector: 'button', match: ['name', 'text']});
  kind('input', {selector: '[data-u2$="-input"]', match: ['name', 'label'], labelSelector: '.u2-input-label'});
  kind('dialog', {selector: '[data-u2="dialog"]', match: ['name', 'title'], labelSelector: '.u2-dialog-title-text'});
  kind('toolbar', {selector: '[data-u2="toolbar"]', match: ['name']});
  element('toolbox', {selector: '[data-u2-name="toolbox"]', aliases: ['workbench toolbar']});
  element('save button inside toolbox', {selector: '[data-u2-name="toolboxSave"]'});
  element('msa dialog', {selector: '[data-u2-name="msaDialog"]', parts: {'title bar': '.u2-dialog-title'}});
});

test('a registered whole phrase wins over composition', () => {
  const ref = parseNoun('save button inside toolbox');
  assert.equal(ref.plan.type, 'entry');
  assert.equal(ref.scope, undefined);
  assert.match(describeNoun(ref), /toolboxSave/);
});

test('composition scopes the inner phrase inside the outer one', () => {
  const ref = parseNoun('run button in the toolbox');
  assert.equal(ref.plan.type, 'kind');
  assert.equal((ref.plan as any).kind.name, 'button');
  assert.equal((ref.plan as any).qualifier, 'run');
  assert.equal(ref.scope?.plan.type, 'entry');
  assert.equal((ref.scope?.plan as any).entry.name, 'toolbox');
});

test('generic kind by longest suffix, qualifier kept', () => {
  const ref = parseNoun('sequence column input in MSA dialog');
  assert.equal((ref.plan as any).qualifier, 'sequence column');
  assert.equal((ref.scope?.plan as any).entry.name, 'msa dialog');
});

test('every kind suffix is kept as an alternative, longest first', () => {
  kind('column input', {selector: '[data-u2="column-input"]', match: ['name']});
  const ref = parseNoun('sequence column input');
  assert.equal((ref.plan as any).kind.name, 'column input');
  assert.equal((ref.plan as any).qualifier, 'sequence');
  assert.deepEqual((ref.plan as any).alternatives.map((a: any) => [a.kind.name, a.qualifier]), [['input', 'sequence column']]);
  assert.match(describeNoun(ref), /column input.*or kind "input" qualified "sequence column"/);
});

test('a phrase that is itself a kind name tries the qualified reading first', () => {
  kind('columns input', {selector: '[data-u2="columns-input"]', match: ['name']});
  const ref = parseNoun('columns input');
  assert.equal((ref.plan as any).kind.name, 'input');
  assert.equal((ref.plan as any).qualifier, 'columns');
  assert.equal((ref.plan as any).alternatives[0].kind.name, 'columns input');
  assert.equal((ref.plan as any).alternatives[0].qualifier, '');
});

test('"of" finds the part on whichever reading of the scope defines it', () => {
  kind('form', {selector: '[data-u2="form"]', match: ['name']});
  kind('function form', {selector: '[data-u2="func-form"]', match: ['name'], parts: {'history icon': '[data-u2="ff-history-icon"]'}});
  const ref = parseNoun('history icon of function form');
  assert.equal(ref.plan.type, 'part');
  assert.equal((ref.plan as any).selector, '[data-u2="ff-history-icon"]');
  assert.equal((ref.scope?.plan as any).kind.name, 'function form');
  assert.equal((ref.scope?.plan as any).alternatives[0].kind.name, 'form');
});

test('a kind alone is an unqualified match', () => {
  const ref = parseNoun('the dialog');
  assert.equal(ref.plan.type, 'kind');
  assert.equal((ref.plan as any).qualifier, '');
});

test('aliases resolve like names', () => {
  assert.equal((parseNoun('the workbench toolbar').plan as any).entry.name, 'toolbox');
});

test('ordinals', () => {
  assert.equal(parseNoun('second save button').ordinal, 1);
  assert.equal(parseNoun('last button in toolbox').ordinal, 'last');
  assert.equal(parseNoun('3rd input').ordinal, 2);
});

test('quoted qualifiers keep scope words intact', () => {
  const ref = parseNoun('"Sign in" button in toolbox');
  assert.equal((ref.plan as any).qualifier, 'sign in');
  assert.equal((ref.scope?.plan as any).entry.name, 'toolbox');
});

test('"of" names a registered part', () => {
  const ref = parseNoun('title bar of MSA dialog');
  assert.equal(ref.plan.type, 'part');
  assert.equal((ref.plan as any).selector, '.u2-dialog-title');
});

test('nested scopes chain', () => {
  const ref = parseNoun('ok button in msa dialog in toolbox');
  assert.equal(ref.scope?.plan.type, 'entry');
  assert.equal(ref.scope?.scope?.plan.type, 'entry');
});

test('unknown phrases fail with the kind list', () => {
  assert.throws(() => parseNoun('the frobnicator'), (e: unknown) => e instanceof NounError && /known kind/.test(e.message));
});

test('a context supplies names of its own; platform names keep their global meaning inside it', () => {
  const workbench = context('MSA workbench', {selector: '[data-u2-name="msaWorkbench"]'});
  workbench.element('results', {selector: '[data-u2-name="results"]'});
  const local = parseNoun('results', workbench);
  assert.equal(local.plan.type, 'entry');
  assert.equal((local.plan as any).entry.context, workbench);
  assert.equal(contextFirst(local), true);
  const global = parseNoun('toolbox', workbench);
  assert.equal((global.plan as any).entry.context, undefined);
  assert.equal(contextFirst(global), false);
  assert.equal(contextFirst(parseNoun('toolbar', workbench)), true);
  assert.equal(contextFirst(parseNoun('save button in toolbar', workbench)), false);
  assert.throws(() => parseNoun('results'), NounError);
  assert.match(describeNoun(local), /of msa workbench/);
  assert.match(describeNoun(parseNoun('toolbar', workbench)), /inside msa workbench first/);
  assert.equal((parseNoun('msa workbench', workbench).plan as any).entry, workbench);
});

test('names are unique: a platform name cannot be registered again, not even on a context', () => {
  assert.throws(() => element('toolbox', {selector: '.other'}), /already registered.*register app names on a context/);
  const workbench = context('MSA workbench', {selector: '.w'});
  assert.throws(() => workbench.element('toolbox', {selector: '.bar'}), /platform names are reserved inside "msa workbench" too/);
  assert.throws(() => context('toolbox', {selector: '.x'}), /already registered/);
  workbench.element('results', {selector: '.r'});
  assert.throws(() => workbench.element('results', {selector: '.r2'}), /already registered/);
});
