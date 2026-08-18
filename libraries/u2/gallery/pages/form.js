import {signal, computed, Scope, Component} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {TextInput, TextArea} from '../../src/components/text-input.js';
import {BoolInput} from '../../src/components/bool-input.js';
import {Form} from '../../src/components/form.js';

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
injectOnce('u2-form-css', '../../css/form.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function sampleForm(options) {
  const form = new Form(options);
  const name = new TextInput({label: 'Compound name', value: 'Aspirin'});
  name.addValidator((v) => v.trim().length > 0 ? null : 'Name is required');
  const description = new TextArea({label: 'Description', placeholder: 'Up to 40 characters'});
  description.addValidator((v) => v.length > 40 ? 'At most 40 characters' : null);
  const active = new BoolInput({label: 'Active', value: true});
  return form.addAll([name, description, active]);
}

const DEFAULTS = {'Compound name': 'Aspirin', 'Description': '', 'Active': true};

export async function render(main) {
  main.append(el('h1', null, 'Form'));
  const intro = el('p');
  intro.innerHTML = 'A form lays inputs out — one label column measured across the rows — and ' +
    'aggregates their <code>validity</code> signals; <code>Enter</code> moves to the next editor. ' +
    'It never owns its inputs: each one is adopted by the scope that created it, so an input may ' +
    'outlive the form that showed it. The P6 schema-driven generator builds on this container.';
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

  const saved = signal('nothing saved yet');
  section('Validation and buttons', () => {
    const form = sampleForm();
    form.addButtons((row) => {
      row.append(
        button('Set values', () => form.setValues(DEFAULTS)),
        button('Save', () => {
          saved.value = form.validate() ? JSON.stringify(form.getValues()) : 'form is invalid';
        }, {primary: true})
      );
    });
    return [
      form,
      readout('validity', computed(() => form.validity.value ?? 'valid')),
      readout('saved', saved),
    ];
  });

  section('Default and condensed, side by side', () => {
    const row = divH([
      divV([span('default'), sampleForm()]),
      divV([span('condensed'), sampleForm({condensed: true})]),
    ]);
    row.style.gap = 'var(--dg-space-xxl)';
    return [row];
  });

  section('Wide (two columns)', () => [sampleForm({wide: true})]);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each form above was built inside a Component.build(...) builder, so the builder owns both the ' +
    'form and its inputs: disposing the sections releases every effect and listener, and live ' +
    'scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
