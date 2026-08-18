/* Wizard: gating, navigation, lazy steps and dialog mode on the node DOM shim. Same contract as
   smoke.test.js — every test leaves the live-scope count where it found it. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal, Scope} from '../src/index.js';
import {Wizard} from '../src/components/wizard.js';

function wizard(name, body) {
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

function markers(w) {
  return w.root.querySelectorAll('.u2-wizard-step');
}

function footer(w, text) {
  return w.root.querySelectorAll('.u2-wizard-buttons button').find((b) => b.textContent === text);
}

function content(title, text) {
  const el = document.createElement('div');
  el.textContent = `${title}: ${text}`;
  return el;
}

function steps(gate) {
  return [
    {id: 'one', title: 'One', content: content('One', 'free')},
    {id: 'two', title: 'Two', content: content('Two', 'gated'), canProceed: gate},
    {id: 'three', title: 'Three', content: content('Three', 'last')},
  ];
}

wizard('navigates forward and back, keeping panel state', () => {
  const w = new Wizard({steps: steps(undefined)});
  document.body.append(w.root);

  const panels = w.root.querySelectorAll('.u2-wizard-panel');
  assert.equal(w.currentStep.value, 'one');
  assert.equal(panels[0].style.display, '');
  assert.equal(panels[1].style.display, 'none');
  assert.equal(markers(w)[0].getAttribute('aria-current'), 'step');
  assert.equal(footer(w, 'BACK').style.display, 'none', 'BACK hidden on the first step');

  w.next();
  assert.equal(w.currentStep.value, 'two');
  assert.equal(footer(w, 'BACK').style.display, '');
  assert.equal(markers(w)[0].classList.contains('u2-wizard-step-done'), true);
  assert.equal(markers(w)[0].textContent.startsWith('✓'), true);

  w.next();
  assert.ok(footer(w, 'FINISH'), 'NEXT becomes FINISH on the last step');
  w.back();
  assert.equal(w.currentStep.value, 'two');
  assert.equal(panels[1].textContent, 'Two: gated', 'panels are hidden, never rebuilt');
  w.dispose();
});

wizard('blocks on canProceed and shows the reason', () => {
  const gate = signal('Fill in the form first');
  const w = new Wizard({steps: steps(gate)});
  document.body.append(w.root);

  w.next();
  assert.equal(w.currentStep.value, 'two');
  const next = footer(w, 'NEXT');
  const reason = w.root.querySelector('.u2-wizard-reason');
  assert.equal(next.disabled, true);
  assert.equal(reason.textContent, 'Fill in the form first');

  w.next();
  assert.equal(w.currentStep.value, 'two', 'a blocked step does not advance');

  gate.value = null;
  assert.equal(next.disabled, false);
  assert.equal(reason.textContent, '');
  w.next();
  assert.equal(w.currentStep.value, 'three');
  w.dispose();
});

wizard('re-evaluates a callback gate on every attempt', () => {
  let allowed = false;
  const w = new Wizard({steps: steps(() => allowed ? true : 'Not yet')});
  document.body.append(w.root);
  w.next();

  w.next();
  assert.equal(w.currentStep.value, 'two');
  allowed = true;
  w.next();
  assert.equal(w.currentStep.value, 'three', 'the callback is re-read on the next attempt');
  w.dispose();
});

wizard('builds lazy content once and fires onActivate on every visit', () => {
  let built = 0;
  const activated = [];
  const w = new Wizard({steps: [
    {id: 'one', title: 'One', content: content('One', 'free'), onActivate: (s) => activated.push(s.id)},
    {
      id: 'two',
      title: 'Two',
      onActivate: (s) => activated.push(s.id),
      content: () => {
        built++;
        return content('Two', 'lazy');
      },
    },
  ]});
  document.body.append(w.root);
  assert.equal(built, 0);

  w.next();
  assert.equal(built, 1);
  w.back();
  w.next();
  assert.equal(built, 1, 'lazy content is built once');
  assert.deepEqual(activated, ['one', 'two', 'one', 'two']);
  w.dispose();
});

wizard('goTo only reaches visited steps or the next allowed one', () => {
  const gate = signal('blocked');
  const w = new Wizard({steps: steps(gate)});
  document.body.append(w.root);

  w.goTo('three');
  assert.equal(w.currentStep.value, 'one', 'unvisited steps beyond the next one are ignored');
  w.goTo('two');
  assert.equal(w.currentStep.value, 'two', 'the immediate next step goes through next()');
  w.goTo('three');
  assert.equal(w.currentStep.value, 'two', 'blocked by the gate, like next()');
  w.goTo('one');
  assert.equal(w.currentStep.value, 'one', 'visited steps are always reachable');

  fire(markers(w)[1], 'click');
  assert.equal(w.currentStep.value, 'two', 'clicking a visited marker navigates');
  w.dispose();
});

wizard('keyboard: Enter advances, arrows move between markers', () => {
  let finished = 0;
  const w = new Wizard({steps: steps(undefined), onFinish: () => finished++});
  document.body.append(w.root);

  const panel = w.root.querySelector('.u2-wizard-panel');
  fire(panel, 'keydown', {key: 'Enter'});
  assert.equal(w.currentStep.value, 'two');
  fire(panel, 'keydown', {key: 'Enter', defaultPrevented: true});
  assert.equal(w.currentStep.value, 'two', 'a handled Enter is left alone');

  markers(w)[1].focus();
  fire(markers(w)[1], 'keydown', {key: 'ArrowLeft'});
  assert.equal(document.activeElement, markers(w)[0], 'focus moves without navigating');
  assert.equal(w.currentStep.value, 'two');
  fire(markers(w)[0], 'keydown', {key: 'Enter'});
  assert.equal(w.currentStep.value, 'one');

  w.goTo('two');
  w.next();
  fire(w.root.querySelectorAll('.u2-wizard-panel')[2], 'keydown', {key: 'Enter', ctrlKey: true});
  assert.equal(finished, 1);
  assert.equal(w.completed.value, true);
  w.dispose();
});

wizard('dialog mode: finishes, cancels on close, and disposes clean', () => {
  const outcome = [];
  const scope = new Scope();
  const w = Scope.runWith(scope, () => new Wizard({
    steps: steps(undefined),
    onFinish: () => outcome.push('finish'),
    onCancel: () => outcome.push('cancel'),
  }));

  const dialog = w.openInDialog('Set up', {width: 400});
  assert.equal(dialog.isOpen.value, true);
  assert.equal(dialog.root.contains(w.root), true);
  assert.equal(dialog.root.querySelectorAll('.u2-dialog-buttons button').length, 0);

  dialog.close();
  assert.deepEqual(outcome, ['cancel'], 'any close that is not a finish cancels');

  w.openInDialog('Set up');
  assert.equal(dialog.isOpen.value, true, 'the same dialog reopens');
  footer(w, 'CANCEL').click();
  assert.deepEqual(outcome, ['cancel', 'cancel']);

  w.openInDialog('Set up');
  w.next();
  w.next();
  footer(w, 'FINISH').click();
  assert.deepEqual(outcome, ['cancel', 'cancel', 'finish']);
  assert.equal(dialog.isOpen.value, false, 'FINISH closes the dialog');

  scope.dispose();
});
