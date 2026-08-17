import {signal, computed, Scope, Component} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {DateInput, DateTimeInput} from '../../src/components/date-input.js';

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
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-date-css', '../../css/date.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function iso(value) {
  return value === null ? 'null' : value.toISOString();
}

export async function render(main) {
  main.append(el('h1', null, 'Date & time'));
  const intro = el('p');
  intro.innerHTML = 'A text editor plus a calendar popup, sharing one hand-written machine: ' +
    '<code>DateInput</code> reads and writes <code>yyyy-MM-dd</code> at midnight local, ' +
    '<code>DateTimeInput</code> adds a 24-hour <code>HH:mm</code> row and zeroes seconds. ' +
    'Values are native <code>Date | null</code>. Typing that does not parse stays on screen and ' +
    'marks the input invalid — the value signal keeps the last good date; min/max clamp on ' +
    'commit (blur, Enter), never mid-keystroke. In the calendar: arrows move the day, ' +
    'PageUp/PageDown change the month (Shift the year), Home/End jump to the week edges, ' +
    'Enter or a click selects, Esc closes. v1 is single-date, local time — no ranges, no locales.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Component.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  const now = new Date();
  const midnight = new Date(now.getFullYear(), now.getMonth(), now.getDate());
  const day = signal(midnight);
  section('Date', () => {
    const input = new DateInput({label: 'Start', bind: day,
      tooltipText: 'yyyy-MM-dd; the calendar icon or ArrowDown opens the picker'});
    return [input, readout('start', computed(() => iso(day.value)))];
  });

  const bounded = signal(midnight);
  section('Bounded (today ±1 week)', () => {
    const input = new DateInput({label: 'Visit', bind: bounded,
      min: new Date(now.getFullYear(), now.getMonth(), now.getDate() - 7),
      max: new Date(now.getFullYear(), now.getMonth(), now.getDate() + 7)});
    return [
      input,
      readout('visit', computed(() => iso(bounded.value))),
      span('Days outside the range are disabled in the calendar; a typed date outside it ' +
        'is clamped on commit.', 'u2-gallery-status'),
    ];
  });

  const when = signal(new Date(now.getFullYear(), now.getMonth(), now.getDate(), 9, 30));
  section('Date and time', () => {
    const input = new DateTimeInput({label: 'Runs at', bind: when});
    return [
      input,
      readout('runsAt', computed(() => iso(when.value))),
      span('The footer row edits HH:mm — ArrowUp/ArrowDown step and wrap, two digits advance ' +
        'to the minutes, Enter commits date and time together. Picking another day keeps the ' +
        'time of day.', 'u2-gallery-status'),
    ];
  });

  const shared = signal(new Date(now.getFullYear(), now.getMonth(), 1));
  section('Two editors, one signal', () => [
    new DateInput({label: 'From', bind: shared}),
    new DateInput({label: 'Same date', bind: shared}),
    readout('from', computed(() => iso(shared.value))),
  ]);

  section('Typed input', () => {
    const input = new DateInput({label: 'Type 2026-02-30', value: new Date(2026, 1, 28)});
    return [
      input,
      readout('validity', computed(() => input.validity.value ?? 'valid')),
      readout('value', computed(() => iso(input.value.value))),
    ];
  });

  section('Week starts on Sunday', () => {
    const input = new DateInput({label: 'Date', firstDayOfWeek: 0, value: midnight});
    return [input, span('firstDayOfWeek: 0 — the default is 1 (Monday).', 'u2-gallery-status')];
  });

  section('Inline / compact', () => {
    const from = new DateInput({inline: true, value: midnight});
    const toolbar = divH([span('Since'), from]);
    toolbar.style.gap = 'var(--dg-space-m)';
    toolbar.style.alignItems = 'center';
    toolbar.style.padding = 'var(--dg-space-s)';
    toolbar.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    toolbar.style.borderRadius = 'var(--dg-radius)';
    return [toolbar];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Component.build(...) owner: disposing it closes any open calendar, ' +
    'releases the editors\' effects, their key/click listeners and the readout bindings, and ' +
    'live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
