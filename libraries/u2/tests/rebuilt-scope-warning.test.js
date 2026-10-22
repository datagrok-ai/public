/* An effect or a cleanup registered on a long-lived control from inside each() rows or effect runs
   piles up with every rebuild: with u2 debug tracking on, u2 warns once per call site. Registering
   on a control built in the row, or calling from the control's own build, stays quiet. */

import {after, before, test} from 'node:test';
import assert from 'node:assert/strict';
import {flush, resetDom} from './dom-shim.js';
import {debug} from '../src/core/scope.js';
import {signal, Scope, Control, TextInput, each, div} from '../src/index.js';

before(() => debug.enable({stacks: false, orphanSeconds: 0}));
after(() => debug.disable());

class Editor extends Control {
  constructor(rows, register) {
    super();
    this.runInScope(() => {
      const host = div([]);
      each(host, rows, (key) => key, (key) => {
        const input = new TextInput({label: key});
        register(this, input);
        return input.root;
      });
      this.root.append(host);
    });
  }
}

async function warningsOf(register) {
  const warnings = [];
  const warn = console.warn;
  console.warn = (...a) => warnings.push(a.join(' '));
  const rows = signal(['R1', 'R2']);
  const editor = new Editor(rows, register);
  try {
    rows.value = ['R1', 'R2', 'R3'];
    rows.value = ['R3'];
    rows.value = ['R1', 'R2', 'R3', 'R4'];
  }
  finally {
    console.warn = warn;
    editor.dispose();
    resetDom();
    await flush();
  }
  return warnings;
}

test('an effect on the editor from its rows warns once, naming the rows', async () => {
  const live = Scope.liveCount;
  const warnings = await warningsOf((editor) => editor.effect(() => {}));
  assert.equal(warnings.length, 1, warnings.join('\n'));
  assert.match(warnings[0], /effect\(\) on .* is called from inside a rebuilt each\(\) row/);
  assert.equal(Scope.liveCount, live);
});

test('a cleanup on the editor from its rows warns too', async () => {
  const warnings = await warningsOf((editor) => editor.own(() => {}));
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /own\(\) on /);
});

test('the same effect on the row\'s own input is quiet', async () => {
  const warnings = await warningsOf((editor, input) => input.effect(() => {}));
  assert.deepEqual(warnings, []);
});

test('a control registering on itself from its own effect runs warns: every run adds one', async () => {
  const warnings = [];
  const warn = console.warn;
  console.warn = (...a) => warnings.push(a.join(' '));
  const n = signal(0);
  const c = new Control();
  try {
    c.effect(() => {
      n.value;
      c.own(() => {});
    });
    n.value = 1;
    n.value = 2;
  }
  finally {
    console.warn = warn;
    c.dispose();
  }
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /from inside a rebuilt effect run/);
});
