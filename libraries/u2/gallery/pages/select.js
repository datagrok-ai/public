import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {ChoiceInput, MultiChoiceInput} from '../../src/components/inputs/choice-input.js';

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
injectOnce('u2-choice-css', '../../css/choice.css');

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
  main.append(el('h1', null, 'Choice inputs'));
  const intro = el('p');
  intro.innerHTML = '<code>ChoiceInput</code> is a native <code>&lt;select&gt;</code> skinned with the ' +
    'platform chrome — the same convention <code>ui.input.choice</code> follows, so the picker is the ' +
    'one the OS draws. Items are plain strings or <code>{value, label}</code> pairs. ' +
    '<code>MultiChoiceInput</code> is a compact checkbox list whose value holds the checked items in ' +
    'item order. Both sit on the shared input foundation: label column, validity, inline variant.';
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

  const color = signal('Green');
  section('Single choice', () => {
    const plain = new ChoiceInput({label: 'Color', items: ['Red', 'Green', 'Blue'], bind: color, nullable: false});
    const solvent = new ChoiceInput({label: 'Solvent', items: ['Water', 'Ethanol', 'DMSO'],
      tooltipText: 'Leave empty to clear'});
    const status = new ChoiceInput({
      label: 'Status',
      items: [{value: 'new', label: 'New'}, {value: 'wip', label: 'In progress'}, {value: 'done', label: 'Done'}],
      value: 'wip',
      nullable: false,
    });
    return [
      plain, solvent, status,
      readout('color', color),
      readout('solvent', computed(() => solvent.value.value ?? '(null)')),
      readout('status', computed(() => status.value.value ?? '(null)')),
    ];
  });

  section('Changing the item list', () => {
    const colors = ['Red', 'Green', 'Blue'];
    const shades = ['Green', 'Teal', 'Navy'];
    const pick = new ChoiceInput({label: 'Item', items: colors, value: 'Green'});
    const row = divH([
      button('Colors', () => pick.setItems(colors)),
      button('Shades', () => pick.setItems(shades)),
    ]);
    row.style.gap = 'var(--dg-space-m)';
    return [
      pick, row,
      readout('value', computed(() => pick.value.value ?? '(null)')),
      el('p', 'u2-gallery-status', 'Green belongs to both lists, so it survives a swap; Red or Navy ' +
        'does not, and the value falls back to null.'),
    ];
  });

  const solvents = signal(['Ethanol']);
  section('Multiple choice', () => {
    const multi = new MultiChoiceInput({label: 'Solvents', items: ['Water', 'Ethanol', 'DMSO', 'Acetone'],
      bind: solvents});
    return [
      multi,
      readout('solvents', computed(() => solvents.value.length ? solvents.value.join(', ') : '(none)')),
      el('p', 'u2-gallery-status', 'Checkboxes are natively focusable: Tab through the list, Space toggles. ' +
        'The value is a new array in item order on every change, never a mutated one.'),
    ];
  });

  section('Validation', () => {
    const assay = new ChoiceInput({label: 'Assay', items: ['IC50', 'EC50', 'Ki']});
    assay.addValidator((v) => v === null ? 'Value is required' : null);
    return [assay, readout('validity', computed(() => assay.validity.value ?? 'valid'))];
  });

  section('Inline / compact', () => {
    const scope = new ChoiceInput({inline: true, items: ['All rows', 'Selected', 'Filtered'], nullable: false,
      value: 'All rows'});
    const toolbar = divH([span('Show'), scope]);
    toolbar.style.gap = 'var(--dg-space-m)';
    toolbar.style.alignItems = 'center';
    toolbar.style.padding = 'var(--dg-space-s)';
    toolbar.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    toolbar.style.borderRadius = 'var(--dg-radius)';
    return [toolbar, readout('scope', computed(() => scope.value.value ?? '(null)'))];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section was built inside a Control.build(...) builder, so its inputs — their effects and ' +
    'change listeners — are owned by it: disposing the sections drops live scopes back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
