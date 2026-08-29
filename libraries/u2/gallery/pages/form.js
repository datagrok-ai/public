import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, divV, span, button} from '../../src/core/elements.js';
import {TextInput, TextArea} from '../../src/components/inputs/text-input.js';
import {BoolInput} from '../../src/components/inputs/bool-input.js';
import {NumberInput} from '../../src/components/inputs/number-input.js';
import {ChoiceInput} from '../../src/components/inputs/choice-input.js';
import {Form} from '../../src/components/forms/form.js';
import {Section} from '../../src/components/containers/section.js';
import {Splitter} from '../../src/components/containers/splitter.js';

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
injectOnce('u2-number-css', '../../css/number.css');
injectOnce('u2-choice-css', '../../css/choice.css');
injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-section-css', '../../css/section.css');
injectOnce('u2-splitter-css', '../../css/splitter.css');

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

/* text + number + choice + bool: the row set that shows the tall exceptions (bool stays inline,
   everything else wraps its caption above) */
function layoutForm(options) {
  return new Form(options).addAll([
    new TextInput({label: 'Compound name', value: 'Aspirin'}),
    new NumberInput({label: 'Dose level', value: 50, min: 0, max: 100, postfix: 'mg'}),
    new ChoiceInput({label: 'Country', items: ['USA', 'France', 'Germany'], value: 'France'}),
    new BoolInput({label: 'Active', value: true}),
  ]);
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
    const component = Control.build(builder);
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

  section('Layouts: normal, wide, tall, side by side', () => {
    const row = divH([
      divV([span('normal'), layoutForm({layout: 'normal'})]),
      divV([span('wide'), layoutForm({layout: 'wide'})]),
      divV([span('tall'), layoutForm({layout: 'tall'})]),
    ]);
    row.style.gap = 'var(--dg-space-xxl)';
    row.style.alignItems = 'flex-start';
    for (const col of row.children)
      col.style.flex = '1';
    return [row];
  });

  section('Auto layout (drag the sash)', () => {
    const auto = layoutForm();
    const left = divV([auto]);
    left.style.padding = 'var(--dg-space-m)';
    const right = divV([span('Auto picks wide or tall by fit: at rest the editors fill the row ' +
      '(the wide layout); squeeze the left panel until the label column leaves less room than ' +
      '--dg-form-min-editor-width and the captions move above, back to wide 10px later.')]);
    right.style.padding = 'var(--dg-space-m)';
    const splitter = new Splitter([left, right], {direction: 'horizontal', sizes: [0.55, 0.45]});
    splitter.root.style.height = '220px';
    return [splitter];
  });

  section('Caption alignment', () => {
    const align = signal('right');
    const form = layoutForm({captionAlign: align});
    return [
      button('Toggle captionAlign', () => align.value = align.value === 'right' ? 'left' : 'right'),
      form,
      readout('captionAlign', align),
    ];
  });

  section('Sections in a form', () => {
    const form = new Form({layout: 'normal'});
    const dosing = new Section({title: 'Dosing'});
    form.addElement(dosing.root);
    form.add(new NumberInput({label: 'Dose level', value: 50, postfix: 'mg'}), dosing.body);
    const study = new Section({title: 'Study'});
    form.addElement(study.root);
    form.add(new ChoiceInput({label: 'Country', items: ['USA', 'France'], value: 'USA'}), study.body);
    form.add(new TextInput({label: 'City', value: 'Boston'}), study.body);
    return [
      span('Hover a header: the chevron appears left of the caption; click collapses the section.'),
      form,
    ];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each form above was built inside a Control.build(...) builder, so the builder owns both the ' +
    'form and its inputs: disposing the sections releases every effect and listener, and live ' +
    'scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
