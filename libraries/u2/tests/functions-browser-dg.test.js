/* The dg layer of the functions browser, stub-style: the `u2-functions-browser` meta registered
   by registerPlatformComponents builds a live control off the doubles' registry via meta.create(),
   and the factory's platform mapping — lowercase tags, tags-plus-role roles, the Dart signature
   string, the scalarOnly / ignorePackages predicates, punctuation-blind ordering — lands in the
   rendered rows. `DG` and `grok` come from tests/dg-stub.mjs. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {Scope} from '../src/core/scope.js';
import {Registry} from '../src/spec/registry.js';
import {Func, Property} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const {registerPlatformComponents} = await import('../src/dg/viewers/registrations.js');
const {functionsBrowser} = await import('../src/dg/entities/functions-browser.js');
const grok = await import('datagrok-api/grok');

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    Func.registry = [
      new Func('MolWeight', {friendlyName: 'Molecular Weight', tags: ['Chem'],
        description: 'Molecular weight of a molecule',
        inputs: [new Property('smiles', 'string')], outputs: [new Property('mw', 'double')]}),
      new Func('Sketch', {namespace: 'Chem', tags: ['chem'],
        outputs: [new Property('view', 'view')]}),
      new Func('_autoPowerGrid', {tags: []}),
      new Func('Sin', {tags: ['Math'], options: {role: 'panel'},
        inputs: [new Property('x', 'double')], outputs: [new Property('y', 'double')]}),
      new Func('HiddenFn', {tags: ['internal']}),
    ];
    try {
      await body();
    } finally {
      Func.registry = [];
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

const rowNames = (fb) => fb.root.querySelectorAll('.u2-list-row[data-u2-func]')
  .map((r) => r.getAttribute('data-u2-func'));

ui('registerPlatformComponents carries the u2-functions-browser meta, and create() builds it', () => {
  const reg = new Registry();
  registerPlatformComponents(reg);
  const meta = reg.get('u2-functions-browser');
  assert.equal(meta.category, 'Data');
  assert.deepEqual(meta.events, ['change', 'activate']);

  const fb = meta.create({showSignature: false, ignorePackages: ['Chem']});
  document.body.append(fb.root);
  assert.equal(fb.root.dataset.u2, 'functions-browser');
  assert.equal(fb.root.classList.contains('u2-fb-no-signature'), true);
  assert.notEqual(fb.root.querySelector('[data-u2="fb-roles"]'), null,
    'the stub functionRoles feed the roles pane');
  assert.equal(fb.root.querySelector('[data-u2="fb-status"]').textContent, '3 of 3',
    'Chem:Sketch is package-ignored, HiddenFn is internal');
  fb.dispose();
});

ui('the factory maps the platform registry: lowercase tags, roles, signature, sort, shell.o', () => {
  const fb = functionsBrowser();
  document.body.append(fb.root);
  fb.root.querySelector('.u2-list').clientHeight = 220;
  fb.query.value = 'x';
  fb.query.value = '';

  const rows = fb.root.querySelectorAll('.u2-list-row[data-u2-func]');
  assert.deepEqual(rows.map((r) => r.getAttribute('data-u2-func')),
    ['_autoPowerGrid', 'MolWeight', 'Sin', 'Chem:Sketch'],
    'ordered by label with punctuation ignored — the underscore entry files under "a"');

  const mw = rows.find((r) => r.getAttribute('data-u2-func') === 'MolWeight');
  assert.equal(mw.querySelector('.u2-fb-label').textContent, 'Molecular Weight');
  assert.equal(mw.querySelector('.u2-fb-sig').textContent, '(smiles) : double');
  assert.equal(mw.title, 'Molecular weight of a molecule');
  const sketch = rows.find((r) => r.getAttribute('data-u2-func') === 'Chem:Sketch');
  assert.equal(sketch.querySelector('.u2-fb-sig').textContent, '() : view',
    'no inputs renders the Dart () form');

  assert.deepEqual([...fb.checkedTags.peek()], []);
  fb.query.value = '#chem';
  assert.deepEqual(rowNames(fb).sort(), ['Chem:Sketch', 'MolWeight'],
    'the Chem tag was lowercased at the boundary');
  fb.query.value = '@panel';
  assert.deepEqual(rowNames(fb), ['Sin'], 'the role option joins the roles');
  fb.query.value = 'molecule';
  assert.deepEqual(rowNames(fb), ['MolWeight'],
    'a plain term reaches the description ("molecule" is not in the name or label)');

  fb.query.value = '';
  fire(fb.root.querySelector('.u2-list-row[data-index="1"]'), 'click');
  assert.equal(grok.shell.o, fb.selected.peek().data, 'selection lands on the shell');
  fb.dispose();
});

ui('scalarOnly keeps one scalar (or dynamic) output; setCurrentObject false leaves the shell alone',
  () => {
    const before = grok.shell.o;
    const fb = functionsBrowser({scalarOnly: true, setCurrentObject: false});
    document.body.append(fb.root);
    fb.root.querySelector('.u2-list').clientHeight = 220;
    fb.query.value = 'x';
    fb.query.value = '';
    assert.deepEqual(rowNames(fb).sort(), ['MolWeight', 'Sin'],
      'the view-typed output and the no-output junk are out');
    fire(fb.root.querySelector('.u2-list-row[data-index="0"]'), 'click');
    assert.equal(grok.shell.o, before);
    fb.dispose();
  });
