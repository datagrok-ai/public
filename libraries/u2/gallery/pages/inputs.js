import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button, bigButton, link} from '../../src/core/elements.js';
import {TextInput, TextArea} from '../../src/components/inputs/text-input.js';
import {BoolInput} from '../../src/components/inputs/bool-input.js';

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

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

export async function render(main) {
  main.append(el('h1', null, 'Inputs & primitives'));
  const intro = el('p');
  intro.innerHTML = 'Fluent element factories (<code>div/divH/span/button/link…</code>) take ' +
    '<code>T | Signal&lt;T&gt;</code> in every text and child slot: a signal argument binds through the ' +
    'ambient <code>Scope</code>, so live values are released with their owner. ' +
    'Inputs share one foundation — label column, underline editor, validity — and every one of them ' +
    'has a compact inline variant.';
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

  const name = signal('world');
  section('Fluent primitives', () => {
    const row = divH([
      button('Say world', () => name.value = 'world'),
      button('Say u2', () => name.value = 'u2'),
      bigButton('Reset', () => name.value = ''),
      link('Callback link', () => name.value = 'a link'),
      link('External link', 'https://datagrok.ai'),
    ]);
    row.style.gap = 'var(--dg-space-m)';
    return [
      readout('greeting', computed(() => `Hello, ${name.value || '…'}!`)),
      row,
    ];
  });

  const query = signal('');
  section('Text', () => {
    const text = new TextInput({label: 'Name', value: 'Aspirin', tooltipText: 'Compound name'});
    const search = new TextInput({label: 'Search', search: true, bind: query,
      placeholder: 'Filter compounds…',
      tooltipText: 'The clear icon shows only while the pointer is over the input, as the platform\'s does'});
    const notes = new TextArea({label: 'Notes', placeholder: 'Grows as you type…', autoGrow: true});
    return [text, search, notes, readout('search', computed(() => query.value || '(empty)'))];
  });

  const notify = signal(true);
  section('Bool', () => [
    new BoolInput({label: 'Enabled', value: true}),
    new BoolInput({label: 'Notifications', switch: true, bind: notify, tooltipText: 'Toggle with Space or Enter'}),
    readout('notifications', computed(() => String(notify.value))),
  ]);

  section('Validation', () => {
    const code = new TextInput({label: 'Code', placeholder: 'up to 10 characters'});
    code.addValidator((v) => v.trim().length > 0 ? null : 'Value is required');
    code.addValidator((v) => v.length > 10 ? 'At most 10 characters' : null);
    return [code, readout('validity', computed(() => code.validity.value ?? 'valid'))];
  });

  section('Inline / compact', () => {
    const filter = new TextInput({inline: true, search: true, placeholder: 'Filter rows…'});
    const selected = new BoolInput({inline: true, switch: true, value: true});
    const enabled = signal(true);
    const toolbar = divH([
      filter,
      selected,
      span('Selected only'),
      button('Toggle enabled', () => {
        enabled.value = !enabled.value;
        filter.enabled = enabled.value;
        selected.enabled = enabled.value;
      }),
    ]);
    toolbar.style.gap = 'var(--dg-space-m)';
    toolbar.style.padding = 'var(--dg-space-s)';
    toolbar.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    toolbar.style.borderRadius = 'var(--dg-radius)';
    return [toolbar, readout('editors enabled', computed(() => String(enabled.value)))];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every input above was created inside a Control.build(...) builder, so it is owned by that ' +
    'component: disposing the sections releases their effects and listeners too, and live scopes ' +
    'drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
