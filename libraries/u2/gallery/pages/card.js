import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {AsyncValue} from '../../src/core/async-value.js';
import {Card} from '../../src/components/containers/card.js';
import {StatCard} from '../../src/components/display/stat-card.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-async-css', '../../css/async.css');
injectOnce('u2-card-css', '../../css/card.css');

const MEDIA = 'data:image/svg+xml;utf8,' + encodeURIComponent(
  '<svg xmlns="http://www.w3.org/2000/svg" width="320" height="120">' +
  '<rect width="320" height="120" fill="#DCEAF8"/>' +
  '<text x="160" y="66" text-anchor="middle" font-family="Roboto,sans-serif" ' +
  'font-size="16" fill="#40607F">C9H8O4</text></svg>');

function el(name, cls, text) {
  const e = document.createElement(name);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function row(children) {
  const r = divH(children);
  r.style.gap = 'var(--dg-space-l)';
  r.style.alignItems = 'flex-start';
  r.style.flexWrap = 'wrap';
  return r;
}

function sized(card, width = '260px') {
  card.root.style.width = width;
  return card;
}

export async function render(main) {
  main.append(el('h1', null, 'Card & stat card'));
  const intro = el('p');
  intro.innerHTML = '<code>Card</code> is the generic surface: header (title, subtitle, leading ' +
    'icon, trailing actions), media, body and footer, each rendered only when given. ' +
    '<code>clickable</code> makes the whole card a button; <code>selectable</code> turns a click ' +
    'into a selection toggle, and a <code>Signal</code> handed as <code>selected</code> is ' +
    'adopted — the card and its owner share one signal. <code>StatCard</code> reuses the same ' +
    'surface for a KPI: a formatted value, its label, an optional signed delta, and — with a ' +
    'source — the skeleton and error states of the <code>AsyncValue</code> contract.';
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

  section('Sections', () => row([
    sized(new Card({
      title: 'Aspirin', subtitle: 'NSAID', icon: 'capsules',
      children: [el('div', null, 'Acetylsalicylic acid, 180.16 g/mol.')],
      footer: [el('div', 'u2-gallery-status', 'ChEMBL25')],
    })),
    sized(new Card({title: 'Header only', subtitle: 'no body content'})),
    sized(new Card({children: [el('div', null, 'Body only — no header, no footer.')]})),
  ]));

  section('Media', () => row([
    sized(new Card({
      title: 'Aspirin', subtitle: 'C9H8O4', media: MEDIA,
      children: [el('div', null, 'The media band takes an image URL, or a ready element.')],
    })),
  ]));

  const clicks = signal(0);
  section('Clickable', () => [
    row([
      sized(new Card({
        title: 'Open the compound', subtitle: 'the whole card is a button', icon: 'external-link',
        children: [el('div', null, 'Click it, or focus it and press Enter or Space.')],
        onClick: () => clicks.value++,
      })),
    ]),
    readout('clicks', computed(() => String(clicks.value))),
  ]);

  const picked = ['Aspirin', 'Caffeine', 'Dopamine'].map(() => signal(false));
  section('Selectable', () => [
    row(['Aspirin', 'Caffeine', 'Dopamine'].map((name, i) => sized(new Card({
      title: name, icon: 'flask', selectable: true, selected: picked[i],
      children: [el('div', null, 'Click to toggle the selection ring.')],
    }), '200px'))),
    row([
      button('Select all', () => picked.forEach((p) => p.value = true)),
      button('Clear', () => picked.forEach((p) => p.value = false)),
    ]),
    readout('selected', computed(() => ['Aspirin', 'Caffeine', 'Dopamine']
      .filter((_, i) => picked[i].value).join(', ') || 'none')),
  ]);

  const revenue = signal(1_240_000);
  const change = signal(0.12);
  section('Stat cards', () => [
    row([
      sized(new StatCard({label: 'Revenue', value: '1.2M', delta: 0.12, icon: 'dollar-sign'}), '200px'),
      sized(new StatCard({
        label: 'Revenue (live)', value: revenue, delta: change, icon: 'chart-line',
        format: (x) => `${(x / 1e6).toFixed(2)}M`,
      }), '200px'),
      sized(new StatCard({label: 'Error rate', value: '0.8%', delta: -0.05, deltaInverted: true,
        icon: 'exclamation-triangle'}), '200px'),
      sized(new StatCard({label: 'Queries today'}), '200px'),
    ]),
    row([
      button('+120k', () => {
        revenue.value += 120_000;
        change.value = 0.12;
      }),
      button('−120k', () => {
        revenue.value -= 120_000;
        change.value = -0.09;
      }),
      button('No change', () => change.value = 0),
      button('Unknown', () => change.value = undefined),
    ]),
    readout('delta', computed(() => change.value === undefined ? 'undefined (the line hides)' :
      `${change.value} (${change.value > 0 ? 'success' : change.value < 0 ? 'failure' : 'neutral'} color)`)),
  ]);

  section('Stat card over a source', () => {
    let outcome = 'ok';
    const source = new AsyncValue(() => new Promise((resolve, reject) =>
      setTimeout(() => outcome === 'ok' ? resolve(48213) : reject(new Error('connection refused')), 900)),
    {auto: false});
    const card = sized(new StatCard({label: 'Queries today', source, icon: 'database'}), '200px');
    const run = (next) => () => {
      outcome = next;
      source.refresh();
    };
    return [
      row([card]),
      row([button('Simulate', run('ok')), button('Fail', run('fail'))]),
      el('div', 'u2-gallery-status',
        'Idle and loading render the skeleton, an error renders its message, ready renders the ' +
        'formatted value — none of it hand-rolled by the card.'),
    ];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every section is a Control.build(...) builder, so the cards constructed inside it — their ' +
    'signal bindings, click listeners and the source above — are owned by it: disposing the ' +
    'sections drops live scopes back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
