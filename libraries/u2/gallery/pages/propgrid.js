import {computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {PropertyGrid} from '../../src/components/property-grid.js';

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
injectOnce('u2-number-css', '../../css/number.css');
injectOnce('u2-property-grid-css', '../../css/property-grid.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

const PROPS = [
  {name: 'Name', type: 'string', description: 'Viewer title shown in the tab header'},
  {name: 'Type', type: 'string', readonly: true, description: 'Viewer type — set at creation'},

  {name: 'Color scheme', type: 'choice', category: 'Appearance', choices: ['Categorical', 'Linear', 'Diverging'],
    description: 'Palette used for the marker color'},
  {name: 'Font size', type: 'int', category: 'Appearance', min: 6, max: 72, description: 'Axis label size, in pixels'},
  {name: 'Opacity', type: 'double', category: 'Appearance', min: 0, max: 1,
    description: 'Marker opacity, clamped to 0…1 on commit'},
  {name: 'Show legend', type: 'bool', category: 'Appearance'},

  {name: 'Width', type: 'int', category: 'Layout', min: 100, max: 4000},
  {name: 'Height', type: 'int', category: 'Layout', min: 100, max: 4000},
  {name: 'Auto layout', type: 'bool', category: 'Layout', description: 'Re-fit the viewer when the window resizes'},

  {name: 'Table', type: 'string', category: 'Data', readonly: true},
  {name: 'Row count', type: 'int', category: 'Data', readonly: true},
  {name: 'Filter', type: 'string', category: 'Data', description: 'Expression applied before rendering'},
];

const VALUES = {
  'Name': 'Scatter plot',
  'Type': 'Scatter plot',
  'Color scheme': 'Categorical',
  'Font size': 12,
  'Opacity': 0.8,
  'Show legend': true,
  'Width': 640,
  'Height': 400,
  'Auto layout': false,
  'Table': 'demog',
  'Row count': 5000,
  'Filter': '',
};

const PRESET = {
  'Name': 'Preset from code',
  'Color scheme': 'Diverging',
  'Font size': 16,
  'Opacity': 0.35,
  'Show legend': false,
  'Width': 1024,
  'Height': 768,
  'Auto layout': true,
  'Row count': 12000,
  'Filter': 'AGE > 40',
};

export async function render(main) {
  main.append(el('h1', null, 'Property grid'));
  const intro = el('p');
  intro.innerHTML = 'Two columns — a name cell (its <code>title</code> carries the description) and an ' +
    'inline editor picked from the property type. Rows sit under collapsible category headers; ' +
    'uncategorized ones come first. Every edit replaces the <code>values</code> record and reports the ' +
    'change through <code>onChanged</code>; writing <code>values</code> from code pushes into all editors ' +
    'without echoing back. Each <code>setProperties</code> builds its rows under a fresh child scope, so ' +
    'no bridge outlives the row it belonged to — while category collapse state survives by name.';
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

  section('Viewer properties', () => {
    const grid = new PropertyGrid();
    grid.setProperties(PROPS, VALUES);
    grid.root.style.width = '380px';
    grid.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    grid.root.style.borderRadius = 'var(--dg-radius)';

    const json = span(computed(() => JSON.stringify(grid.values.value, null, 1)), 'u2-gallery-status');
    json.style.whiteSpace = 'pre-wrap';
    json.style.fontFamily = 'var(--dg-font-family-mono)';

    const lastEdit = computed(() => {
      const change = grid.onChanged.value;
      return change === null ? 'nothing edited yet' : `${change.name} = ${JSON.stringify(change.value)}`;
    });

    const buttons = divH([
      button('Set programmatically', () => grid.values.value = {...grid.values.peek(), ...PRESET}),
      button('Reset', () => grid.values.value = {...VALUES}),
      button('Re-apply properties', () => grid.setProperties(PROPS, grid.values.peek())),
      button('Add a category', () => grid.setProperties(
        PROPS.concat([{name: 'Grid lines', type: 'bool', category: 'Axes'},
          {name: 'Ticks', type: 'int', category: 'Axes', min: 0, max: 20}]),
        grid.values.peek())),
    ]);
    buttons.style.gap = 'var(--dg-space-m)';
    buttons.style.margin = 'var(--dg-space-m) 0';

    return [grid, buttons, readout('last edit', lastEdit), json];
  });

  main.append(el('h2', null, 'Collapsing'));
  main.append(el('p', 'u2-gallery-status',
    'Click a category header (or focus it and press Enter or Space) to collapse it. Collapse a ' +
    'category, then press "Re-apply properties" or "Add a category": the rows are rebuilt from scratch, ' +
    'yet the categories you collapsed stay collapsed because the state is kept by name.'));

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'The section owns the grid, the grid owns one child scope per setProperties call, and that scope ' +
    'owns every editor and both directions of its value bridge. Disposing the section releases all of ' +
    'them, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
