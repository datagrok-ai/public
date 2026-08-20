import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {DynamicInput} from '../../src/components/dynamic-input.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-entity-css', '../../css/entity.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

/* Stands in for the platform's handler card: the gallery is platform-free, and this is exactly
   the seam u2/dg's metaInput fills with ObjectHandler.renderCard. */
function mockCard(user) {
  const card = el('div', 'u2-entity-card');
  card.append(el('div', null, user.name), el('div', 'u2-gallery-status', user.email));
  return card;
}

export async function render(main) {
  main.append(el('h1', null, 'Dynamic input'));
  const intro = el('p');
  intro.innerHTML = 'An input that <em>shows</em> its value instead of editing it — the port of ' +
    'Dart\'s <code>DynamicInput</code>. Every write swaps the rendered element for a fresh one ' +
    'from <code>render</code>; there is no user edit to report, so the value only ever changes ' +
    'from code. In the platform, <code>u2/dg</code>\'s <code>metaInput()</code> is this input with ' +
    'the object handler\'s card as its renderer (Dart\'s <code>MetaInput</code>); here the card is ' +
    'mocked, since the gallery loads no platform.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Control.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const users = [
    {name: 'Ada Lovelace', email: 'ada@datagrok.ai'},
    {name: 'Grace Hopper', email: 'grace@datagrok.ai'},
  ];
  const current = signal(users[0]);
  section('Rendered object', () => {
    const input = new DynamicInput({label: 'Owner', bind: current, render: (user) =>
      user === null ? null : mockCard(user)});
    const swap = button('Next owner', () =>
      current.value = users[(users.indexOf(current.peek()) + 1) % users.length]);
    const clear = button('Clear', () => current.value = null);
    return [input, divH([swap, clear]),
      readout('owner', computed(() => current.value?.name ?? 'null'))];
  });

  section('Fixed size', () => {
    const input = new DynamicInput({label: 'Preview', size: {width: 240, height: 60},
      value: users[1], render: (user) => mockCard(user)});
    return [input, span('The editor box is sized, the renderer fills it — Dart\'s ' +
      'DynamicInput.size.', 'u2-gallery-status')];
  });

  section('Empty and invalid', () => {
    const input = new DynamicInput({label: 'Owner', render: (user) =>
      user === null ? null : mockCard(user)});
    input.addValidator((value) => value === null ? 'Nothing selected' : null);
    return [divV([input, button('Pick Ada', () => input.value.value = users[0])])];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Disposing a section releases the render effect and detaches whatever the renderer last ' +
    'produced; live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
