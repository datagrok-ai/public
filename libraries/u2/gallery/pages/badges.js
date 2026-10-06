import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {badge, countBadge, dot, tag} from '../../src/components/display/badge.js';
import {ProgressBar} from '../../src/components/display/progress.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-badge-css', '../../css/badge.css');
injectOnce('u2-progress-css', '../../css/progress.css');

const VARIANTS = ['default', 'accent', 'success', 'warning', 'error'];
const STAGES = ['draft', 'review', 'published', 'archived'];

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
  r.style.gap = 'var(--dg-space-m)';
  r.style.alignItems = 'center';
  r.style.flexWrap = 'wrap';
  return r;
}

export async function render(main) {
  main.append(el('h1', null, 'Badges & progress'));
  const intro = el('p');
  intro.innerHTML = 'Status pills, unread counters, status dots and removable tags, plus the ' +
    'progress bar. Text and count slots take <code>T | Signal&lt;T&gt;</code>: a signal binds ' +
    'through the ambient <code>Scope</code>, so the counter follows its source and stops when the ' +
    'owner is disposed. <code>ProgressBar</code> is a component — its value signal is clamped to ' +
    '0..1 and drives the fill, and the indeterminate variant is switchable at runtime.';
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

  const stage = signal('draft');
  section('Badges', () => [
    row(VARIANTS.map((variant) => badge(variant, {variant}))),
    row([
      span('Live text:'),
      badge(computed(() => stage.value.toUpperCase()), {variant: 'accent'}),
      button('Next stage', () => stage.value = STAGES[(STAGES.indexOf(stage.value) + 1) % STAGES.length]),
    ]),
  ]);

  const unread = signal(0);
  section('Count badges', () => {
    const bump = (delta) => () => unread.value = Math.max(0, unread.value + delta);
    return [
      row([
        span('Inbox'), countBadge(unread),
        span('Alerts'), countBadge(unread, {max: 9}),
        span('Static'), countBadge(7),
      ]),
      row([
        button('−1', bump(-1)),
        button('+1', bump(1)),
        button('+10', bump(10)),
        button('100', () => unread.value = 100),
        button('Reset', () => unread.value = 0),
      ]),
      readout('count', computed(() => `${unread.value} (0 hides the badge, >99 reads 99+)`)),
    ];
  });

  section('Status dots', () => row([
    dot(), span('Unknown'),
    dot('accent'), span('Running'),
    dot('success'), span('Passed'),
    dot('warning'), span('Flaky'),
    dot('error'), span('Failed'),
  ]));

  const removed = signal('—');
  section('Tags', () => {
    const chip = (name) => {
      const t = tag(name, {onRemove: () => {
        removed.value = name;
        t.remove();
      }});
      return t;
    };
    return [
      row(['Aspirin', 'Caffeine', 'Dopamine', 'Nicotine'].map(chip).concat([tag('Read-only')])),
      readout('last removed', removed),
    ];
  });

  section('Progress — determinate', () => {
    const bar = new ProgressBar({showPercent: true});
    bar.root.style.width = '360px';
    let timer = 0;
    const simulate = () => {
      bar.value.value = 0;
      clearTimeout(timer);
      const step = () => {
        bar.value.value += 0.05;
        if (bar.value.value < 1)
          timer = setTimeout(step, 120);
      };
      timer = setTimeout(step, 120);
    };
    bar.own(() => clearTimeout(timer));
    return [
      bar,
      row([
        button('+10%', () => bar.value.value += 0.1),
        button('Simulate work', simulate),
        button('Reset', () => {
          clearTimeout(timer);
          bar.value.value = 0;
        }),
      ]),
      readout('value', computed(() => bar.value.value.toFixed(2))),
    ];
  });

  section('Progress — indeterminate', () => {
    const bar = new ProgressBar({indeterminate: true});
    bar.root.style.width = '360px';
    const state = signal('indeterminate');
    return [
      bar,
      row([
        button('Toggle indeterminate', () => {
          bar.indeterminate = !bar.indeterminate;
          state.value = bar.indeterminate ? 'indeterminate' : 'determinate';
          if (!bar.indeterminate)
            bar.value.value = 0.35;
        }),
      ]),
      readout('mode', computed(() => `${state.value} (aria-valuenow is dropped while indeterminate)`)),
    ];
  });

  section('Progress — percent and description', () => {
    const done = signal(12);
    const bar = new ProgressBar({
      showPercent: true,
      description: computed(() => `Uploading ${done.value} of 24 files`),
    });
    bar.root.style.width = '360px';
    bar.effect(() => bar.value.value = done.value / 24);
    const move = (delta) => () => done.value = Math.max(0, Math.min(24, done.value + delta));
    return [
      bar,
      row([button('−6 files', move(-6)), button('+6 files', move(6))]),
      el('div', 'u2-gallery-status', 'The value signal clamps on write, so overshooting stops at 0 and 1.'),
    ];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every section is a Control.build(...) builder, so its bindings, listeners and the simulate ' +
    'timer are owned by it: disposing the sections stops the counters and the running bar, and ' +
    'live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
