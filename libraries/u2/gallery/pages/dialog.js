import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {TextInput} from '../../src/components/inputs/text-input.js';
import {Dialog} from '../../src/components/containers/dialog.js';

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
injectOnce('u2-dialog-css', '../../css/dialog.css');

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
  main.append(el('h1', null, 'Dialog'));
  const intro = el('p');
  intro.innerHTML = 'Dialogs render into the shared <code>Overlay</code> layer: modal ones lay a backdrop ' +
    'just below themselves, and every <code>show()</code> takes the next z offset, so they stack. ' +
    'The title bar drags with pointer capture and keeps the dialog inside the viewport; Esc cancels, ' +
    'Enter confirms (outside textareas and buttons), Tab cycles within the dialog, and focus returns ' +
    'to wherever it was when the dialog opened.';
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

  section('Modal', () => {
    const name = signal('Aspirin');
    const result = signal('(nothing yet)');
    const dialog = Dialog.create('Rename compound')
      .add('Pick a new name — OK and CANCEL both close the dialog and report below.')
      .add(new TextInput({label: 'Name', bind: name}))
      .onOK(() => result.value = `OK: ${name.value}`)
      .onCancel(() => result.value = 'CANCEL');
    return [
      divH([button('Open modal dialog', () => dialog.show({modal: true, width: 360}))]),
      readout('result', result),
      readout('isOpen', computed(() => String(dialog.isOpen.value))),
    ];
  });

  section('Non-modal, positioned, draggable', () => {
    const alias = signal('');
    const applied = signal('(not applied)');
    const dialog = Dialog.create('Properties')
      .add('Non-modal: the page below stays interactive. Drag me by the title bar.')
      .add(new TextInput({label: 'Alias', bind: alias}))
      .addButton('Clear', () => alias.value = '')
      .onOK(() => applied.value = alias.value || '(empty)')
      .onCancel(() => applied.value = '(cancelled)');
    return [
      divH([button('Open at (60, 200)', () => dialog.show({x: 60, y: 200, width: 320}))]),
      readout('applied alias', applied),
      readout('isOpen', computed(() => String(dialog.isOpen.value))),
    ];
  });

  section('Stacked', () => {
    const last = signal('(none)');
    const second = Dialog.create('Second dialog')
      .add('On top of the first one: its backdrop dims the dialog below.')
      .onOK(() => last.value = 'second OK');
    const first = Dialog.create('First dialog')
      .add('Open the second one — each show() takes the next z offset.')
      .addButton('Open second…', () => second.show({modal: true, width: 280}), {primary: true})
      .onOK(() => last.value = 'first OK');
    return [
      divH([button('Open two dialogs', () => {
        first.show({modal: true, width: 380});
        second.show({modal: true, width: 280});
      })]),
      readout('last OK', last),
      readout('first / second open', computed(() => `${first.isOpen.value} / ${second.isOpen.value}`)),
    ];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Every dialog above was created inside a Control.build(...) builder, so the section owns it: ' +
    'disposing closes an open dialog, removes it and its backdrop from the overlay layer, and drops ' +
    'its listeners — live scopes return to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
