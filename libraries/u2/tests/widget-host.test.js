/* host(): a u2 component hosted as a real DG.Widget (DD7). `DG.Widget`, `DG.Property` and
   `grok.functions.register` come from tests/dg-stub.mjs — the platform contract, nothing more;
   everything under test is the shipped module. The compile-time half of the convergence (that
   DG.Widget satisfies WidgetLike) is asserted by tsc in src/dg/widget-host.ts. */

import {test} from 'node:test';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {signal} from '../src/core/signals.js';
import {Scope} from '../src/core/scope.js';
import {Control} from '../src/core/component.js';

register('./dg-stub.mjs', import.meta.url);
const {host, U2Widget} = await import('../src/dg/widget-host.js');
const DG = await import('datagrok-api/dg');
const grok = await import('datagrok-api/grok');

function hosted(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    } finally {
      grok.functions.registrations.length = 0;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

class Gauge extends Control {
  constructor() {
    super();
    this.value = signal(10);
    this.meta = {
      tag: 'u2-gauge',
      description: 'A gauge',
      usage: 'Shows one number.',
      events: ['bumped'],
      props: [{name: 'value', type: 'int'}],
    };
    this.registerFunction({
      name: 'reset',
      description: 'Sets the value back to zero.',
      inputs: [{name: 'to', type: 'int'}],
      apply: (params) => {
        this.value.value = params?.to ?? 0;
        return this.value.peek();
      },
    });
  }

  bump() {
    this.value.value = this.value.peek() + 1;
    this.fireEvent('bumped', this.value.peek());
  }
}

hosted('the host is a platform widget over the component root', () => {
  const g = new Gauge();
  const w = host(g);
  assert.ok(w instanceof DG.Widget);
  assert.ok(w instanceof U2Widget);
  assert.equal(w.root, g.root);
  assert.equal(g.root.getAttribute('data-kill-on-close'), 'true');
  assert.notEqual(w.dart, undefined, 'registered in the platform widget registry');
  w.detach();
});

hosted('detaching the widget disposes the component', () => {
  const g = new Gauge();
  const w = host(g);
  assert.equal(g.scope.isDisposed, false);
  w.detach();
  assert.equal(g.scope.isDisposed, true);
});

hosted('a docked component is detached when its pane is closed', () => {
  const g = new Gauge();
  let onClosed;
  const dockManager = {onClosed: {subscribe: (fn) => {
    onClosed = fn;
    return {unsubscribe() {}};
  }}};
  host(g, dockManager);

  onClosed(document.createElement('div'));
  assert.equal(g.scope.isDisposed, false, 'another pane closing is none of its business');
  onClosed(g.root);
  assert.equal(g.scope.isDisposed, true);
});

hosted('properties are the component\'s own, as platform properties', () => {
  const g = new Gauge();
  const w = host(g);
  const [value] = w.getProperties();
  assert.ok(value instanceof DG.Property);
  assert.equal(value.name, 'value');
  assert.equal(value.propertyType, 'int');
  assert.equal(value.get(w), 10);
  g.bump();
  assert.equal(value.get(w), 11);
  value.set(w, 3);
  assert.equal(g.value.peek(), 3);
  assert.equal(w.getProperties()[0], value, 'minted once');
  w.detach();
});

hosted('functions are minted as platform functions, once, and reach the component', async () => {
  const g = new Gauge();
  const w = host(g);
  const [reset] = w.getFunctions();
  assert.equal(grok.functions.registrations.length, 1);

  const registration = grok.functions.registrations[0];
  assert.equal(registration.signature, 'object u2_gauge_reset(int to)');
  assert.equal(registration.namespace, 'U2');
  assert.equal(registration.isAsync, true);
  assert.equal(registration.options.description, 'Sets the value back to zero.');

  assert.equal(await registration.run(7), 7, 'params arrive positionally, in declared order');
  assert.equal(g.value.peek(), 7);
  assert.equal(w.getFunctions()[0], reset, 'minted once');
  w.detach();
});

hosted('a component with no functions registers nothing', () => {
  const plain = new Control();
  const w = host(plain);
  assert.deepEqual(w.getFunctions(), []);
  assert.equal(grok.functions.registrations.length, 0);
  w.detach();
});

hosted('events reach the platform as an observable', () => {
  const g = new Gauge();
  const w = host(g);
  const seen = [];
  const sub = w.onEvent('bumped').subscribe((x) => seen.push(x));
  const clicks = w.onEvent('click').subscribe((e) => seen.push(e.type));

  g.bump();
  fire(g.root, 'click');
  assert.deepEqual(seen, [11, 'click']);

  sub.unsubscribe();
  clicks.unsubscribe();
  g.bump();
  fire(g.root, 'click');
  assert.deepEqual(seen, [11, 'click']);
  w.detach();
});

hosted('the AI briefing and the status snapshot are the component\'s', () => {
  const g = new Gauge();
  const w = host(g);
  assert.equal(w.aiDescription, 'Shows one number.');
  w.aiDescription = 'Set by the host';
  assert.equal(g.aiDescription, 'Set by the host');

  const status = w.getWidgetStatus();
  assert.deepEqual(status.events, [{name: 'bumped', eventName: 'bumped', description: ''}]);
  assert.equal(status.description, 'A gauge');
  w.detach();
});
