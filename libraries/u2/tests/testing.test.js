import {test} from 'node:test';
import assert from 'node:assert/strict';
import './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope, TextInput, ChoiceInput, button, iconButton, div, signal} from '../src/index.js';
import {u2test, mount, byName, byKind, state, expectOperable, user} from '../src/testing/index.js';

const NO_CASCADE = typeof CSSStyleSheet !== 'function' && 'needs a DOM with the CSS cascade';

/** A page rule cascading over u2's installed sheet, as a package stylesheet would. */
function pageRule(css) {
  const style = document.createElement('style');
  style.textContent = css;
  document.head.append(style);
  return style;
}

function campaignName(onEdit = () => {}) {
  const input = mount(() => new TextInput({label: 'Campaign Name', value: 'KHL-3', readOnly: true}));
  const edit = iconButton('pencil', onEdit, {tooltip: 'Edit Campaign Name'});
  input.addOptions(edit);
  return {input, edit};
}

const OLD_RAIL = '.u2-input-readonly .u2-input-options { pointer-events: none; }';

test('user.click: a caller option on a read-only input is pressed under today\'s rail rule', u2test(() => {
  let edits = 0;
  const {edit} = campaignName(() => edits++);
  user.click(byName(document, 'Edit Campaign Name'));
  assert.equal(edits, 1);
  assert.equal(document.activeElement, edit, 'the press focused the button, as a browser does');
}));

test('user.click: under the old read-only rail rule it is refused with the probe\'s finding, the handler not run',
  {skip: NO_CASCADE}, u2test(() => {
    let edits = 0;
    campaignName(() => edits++);
    const rule = pageRule(OLD_RAIL);
    try {
      assert.throws(() => user.click(byName(document, 'Edit Campaign Name')), {message: 'user.click: button ' +
        '"Edit Campaign Name" takes no pointer: pointer-events: none on div.u2-input-options'});
    }
    finally {
      rule.remove();
    }
    assert.equal(edits, 0, 'a synthetic click would have run it; a user cannot');
  }));

test('user.click: a covered element is refused where the page has layout', u2test(() => {
  let edits = 0;
  const {input, edit} = campaignName(() => edits++);
  layout(edit, {rect: new DOMRect(100, 10, 16, 28)});
  const row = input.root.querySelector('.u2-input-box');
  document.elementFromPoint = (x, y) => x === 108 && y === 24 ? row : null;
  try {
    assert.throws(() => user.click(edit), {message: 'user.click: button "Edit Campaign Name" is covered by ' +
      'div.u2-input-box'});
    document.elementFromPoint = () => edit.firstChild;
    user.click(edit);
  }
  finally {
    delete document.elementFromPoint;
  }
  assert.equal(edits, 1, 'only the uncovered press reached the handler');
}));

test('user.click: off, hidden, inert and detached elements are refused, saying why', u2test(() => {
  let runs = 0;
  const run = () => runs++;
  const why = mount(() => button('Submit', run, {disabled: 'Choose a submit function'}));
  const bare = mount(() => Object.assign(button('Save', run), {disabled: true}));
  const pane = mount(() => div([button('Run', run)], 'pane'));
  const folded = mount(() => div([button('Fold', run)], 'folded'));
  assert.throws(() => user.click(why), {message: 'user.click: button "Submit" is disabled: Choose a submit function'});
  assert.throws(() => user.click(bare), {message: 'user.click: button "Save" is disabled and says no reason'});
  pane.setAttribute('inert', '');
  pane.title = 'Check Compute first';
  assert.throws(() => user.click(byName(pane, 'Run')), {message: 'user.click: button "Run" is in div.pane, which is ' +
    'inert: Check Compute first'});
  folded.hidden = true;
  assert.throws(() => user.click(folded.firstChild), {message: 'user.click: button "Fold" is hidden (in div.folded)'});
  folded.hidden = false;
  folded.style.display = 'none';
  assert.throws(() => user.click(folded.firstChild), {message: 'user.click: button "Fold" is not shown: display: none ' +
    '(in div.folded)'});
  assert.throws(() => user.click(button('Lost', run)), {message: 'user.click: button "Lost" is not on the page'});
  const fenced = mount(() => div([button('Fenced', run)], 'fence'));
  fenced.style.pointerEvents = 'none';
  assert.throws(() => user.click(fenced.firstChild), {message: 'user.click: button "Fenced" takes no pointer: ' +
    'pointer-events: none on div.fence'});
  user.press(fenced.firstChild, 'Enter');
  assert.equal(runs, 1, 'the keyboard still reaches what the pointer cannot');
}));

test('user.type: a read-only field refuses the keys; an editable one takes them as the user\'s', u2test(() => {
  const typed = [];
  const name = mount(() => new TextInput({label: 'Name', value: 'KHL', readOnly: true,
    onInput: (v) => typed.push(v)}));
  assert.throws(() => user.type(byName(document, 'Name'), '-4'), {message: 'user.type: input "Name" is read-only: ' +
    'a user cannot type into it'});
  assert.equal(state(name).value, 'KHL');
  name.readOnly = false;
  user.type(name, '-4');
  assert.equal(state(name).value, 'KHL-4');
  assert.deepEqual(typed, ['KHL-', 'KHL-4'], 'onInput heard every key: the edit counts as the user\'s');
  user.type(name, 'X', {replace: true});
  assert.equal(state(name).value, 'X');
}));

test('user.press: Enter reaches an input\'s onEnter and presses a button', u2test(() => {
  let entered = 0;
  let saved = 0;
  const name = mount(() => new TextInput({label: 'Name', onEnter: () => entered++}));
  const save = mount(() => button('Save', () => saved++));
  user.press(byName(name, 'Name'), 'Enter');
  user.press(save, 'Enter');
  assert.deepEqual([entered, saved], [1, 1]);
  assert.equal(document.activeElement, save);
}));

test('user.pick: an item of a choice input, as the user\'s pick; an item it lacks is refused listing the offer',
  u2test(() => {
    const picked = [];
    const choice = mount(() => new ChoiceInput({label: 'Stage', name: 'stage', items: ['Idea', 'Made'],
      onInput: (v) => picked.push(v)}));
    user.pick(byName(document, 'stage'), 'Made');
    assert.equal(state(choice).value, 'Made');
    assert.deepEqual(picked, ['Made']);
    assert.throws(() => user.pick(choice, 'Tested'), {message: 'user.pick: select "Stage" offers no "Tested"; it ' +
      'offers "", "Idea", "Made"'});
  }));

test('expectOperable: a clean form passes; the old read-only rail rule and a reasonless disable are each listed',
  {skip: NO_CASCADE}, u2test(() => {
    const form = mount(() => div([], 'form'));
    const {input} = campaignName();
    form.append(input.root, button('Save', () => {}, {disabled: 'Name it first'}));
    expectOperable(form);
    form.append(Object.assign(button('Submit', () => {}), {disabled: true}));
    const rule = pageRule(OLD_RAIL);
    try {
      assert.throws(() => expectOperable(form), {message: 'expectOperable: 2 control(s) in div.form cannot be ' +
        'operated as offered:\n  - no-pointer: button "Edit Campaign Name" takes no pointer: pointer-events: none on ' +
        'div.u2-input-options\n  - disabled-without-reason: button "Submit" is disabled and says no reason'});
    }
    finally {
      rule.remove();
    }
    assert.throws(() => expectOperable(div([button('Run', () => {})], 'pane')),
      {message: 'expectOperable: div.pane is not on the page: nothing in it can be operated'});
  }));

test('byName: the automation id, then the name a user reads; none or several say what is there', u2test(() => {
  const {input, edit} = campaignName();
  const status = mount(() => new TextInput({label: 'Status', name: 'status'}));
  const ok = mount(() => button('Submit', () => {}));
  assert.equal(byName(document, 'status'), status.root);
  assert.equal(byName(document, 'edit campaign name'), edit, 'an icon button by its tooltip, case aside');
  assert.equal(byName(document, 'Campaign Name'), input.root.querySelector('input'), 'a field by its label');
  assert.equal(byName(document, 'SUBMIT'), ok, 'a button by its text, as the screen shows it');
  assert.throws(() => byName(document, 'Cancel'), {message: 'byName: nothing in the document is named "Cancel"; ' +
    'named there: status, "Campaign Name", "Edit Campaign Name", "Status", "Submit"'});
  mount(() => button('Submit', () => {}));
  assert.throws(() => byName(document, 'Submit'), {message: 'byName: 2 elements in the document read "Submit": ' +
    'button.u2-btn, button.u2-btn; give the one you mean a name'});
  const hidden = mount(() => button('Cancel', () => {}));
  hidden.hidden = true;
  const shown = mount(() => button('Cancel', () => {}));
  assert.equal(byName(document, 'Cancel'), shown, 'of two that read alike, the one a user sees');
}));

test('byKind and state: control roots by kind; the state a test reads', u2test(() => {
  const {input} = campaignName();
  assert.deepEqual(byKind(document, 'text-input').map((el) => el === input.root), [true]);
  assert.equal(state(input).label, 'Campaign Name');
  assert.equal(state(byName(document, 'Campaign Name')).readOnly, true, 'any element inside the control');
  assert.equal(state(input).value, 'KHL-3');
  assert.throws(() => byKind(document, 'textinput'), {message: 'byKind: no textinput in the document (u2 has no kind ' +
    '"textinput"); kinds there: text-input, icon-button, icon'});
  assert.throws(() => state(document.body), {message: 'state: body is not inside a u2 control'});
}));

test('mount: attached under the test, detached and disposed when it ends; outside a u2test it refuses', async () => {
  let input;
  let bare;
  await u2test(() => {
    input = mount(() => new TextInput({label: 'Name'}));
    bare = mount(() => document.createElement('section'));
    assert.equal(input.root.isConnected, true);
    assert.equal(bare.isConnected, true);
  })();
  assert.equal(input.root.isConnected, false);
  assert.equal(input.scope.isDisposed, true);
  assert.equal(bare.isConnected, false, 'an element no dispose detaches goes too');
  assert.throws(() => mount(() => new TextInput({label: 'Name'})),
    {message: 'mount runs inside a u2test body: test(name, u2test(() => …))'});
});

test('mount: after an await too, what it builds is the test\'s', u2test(async () => {
  await Promise.resolve();
  const late = mount(() => new TextInput({label: 'Late'}));
  assert.equal(late.root.isConnected, true);
}));

test('u2test: a scope the body leaves live fails the test, naming it', async () => {
  let stray;
  await assert.rejects(u2test(() => {
    stray = Scope.root('stray');
  })(), /u2test: the test left 1 scope\(s\), 0 timer\(s\) and 0 global listener\(s\) live once its own scope was disposed: scope stray/);
  stray.dispose();
});

test('u2test: a timer the body leaves pending fails the test', async () => {
  let timer;
  await assert.rejects(u2test(() => {
    const owner = Scope.root('poll');
    owner.build(() => timer = setInterval(() => {}, 1000));
    owner.dispose();
  })(), /the test left 0 scope\(s\), 1 timer\(s\) and 0 global listener\(s\) live .*: interval 1000 ms/);
  clearInterval(timer);
});

test('u2test: a u2 usage error raised meanwhile fails the test', async () => {
  let unowned;
  await assert.rejects(u2test(async () => {
    await Promise.resolve();
    unowned = new TextInput({label: 'Built after an await'});
    unowned.dispose();
  })(), /u2 usage: owner-missing: TextInput was built with no owner/);
  unowned?.dispose();
});

test('u2test: an error that reached Scope.onUnhandledError fails the test', async () => {
  const logged = [];
  const log = Scope.onUnhandledError;
  Scope.onUnhandledError = (e) => logged.push(e);
  try {
    await assert.rejects(u2test(() => {
      const owner = Scope.root('effects');
      owner.effect(() => {
        throw new Error('boom');
      });
      owner.dispose();
    })(), {message: 'u2test: an error reached Scope.onUnhandledError (a handler, an effect or a cleanup threw): boom'});
  }
  finally {
    Scope.onUnhandledError = log;
  }
  assert.equal(logged.length, 1, 'and it still went where unhandled errors go');
});

test('u2test: the body\'s own failure is what the test reports, its leftovers cleaned', async () => {
  const live = Scope.liveCount;
  let input;
  await assert.rejects(u2test(() => {
    input = mount(() => new TextInput({label: 'Name'}));
    assert.equal(signal(1).peek(), 2);
  })(), {name: 'AssertionError'});
  assert.equal(input.root.isConnected, false);
  assert.equal(Scope.liveCount, live);
});
