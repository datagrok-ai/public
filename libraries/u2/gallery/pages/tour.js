import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button, label} from '../../src/core/elements.js';
import {Tour} from '../../src/components/display/tour.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-buttons-css', '../../css/buttons.css');
injectOnce('u2-tour-css', '../../css/tour.css');

function el(name, cls, text) {
  const e = document.createElement(name);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(caption, source, id) {
  const line = divH([span(`${caption} = `), span(source)], 'u2-gallery-status');
  if (id) line.id = id;
  return line;
}

function row(children) {
  const r = divH(children);
  r.style.gap = 'var(--dg-space-m)';
  r.style.alignItems = 'center';
  r.style.flexWrap = 'wrap';
  return r;
}

function named(element, name) {
  element.dataset.u2Name = name;
  return element;
}

export async function render(main) {
  main.append(el('h1', null, 'Tour'));
  const intro = el('p');
  intro.innerHTML = 'A guided walkthrough over UI that already exists: one dimmed layer with ' +
    'SVG-punched spotlight holes, and a popup anchored to the current target through ' +
    '<code>floating-ui</code> — nothing polls, and scrolling or moving the target moves both. ' +
    'The dim layer is click-through, so a step can ask you to actually use the control it points ' +
    'at. Steps are plain data: <code>{target, content, position?, extra?, advanceOn?}</code>, ' +
    'where a string target is a spec name looked up as <code>[data-u2-name]</code> inside ' +
    '<code>root</code>. <b>Esc</b> skips the tour, <b>Enter</b> advances (unless you are typing ' +
    'in a field), and <b>Tab</b> cycles inside the popup.';
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

  const hasRun = signal(false);
  const result = signal('—');
  let tour = null;
  let runButton = null;
  let statusBar = null;

  const app = section('A mock mini-app', () => {
    runButton = named(button('Run', () => hasRun.value = true), 'runButton');
    statusBar = named(span(computed(() => hasRun.value ? 'ran once' : 'idle')), 'statusBar');
    const nameField = named(el('input'), 'nameInput');
    nameField.placeholder = 'Model name';
    return [
      row([named(button('New', () => {}), 'newButton'), named(button('Open', () => {}), 'openButton'), runButton]),
      row([label('Name'), nameField, statusBar]),
    ];
  });

  section('Run the tour', () => [
    row([
      button('Start tour', () => {
        hasRun.value = false;
        result.value = 'running…';
        tour = Tour.run({
          root: app.root,
          steps: [
            {target: 'nameInput', content: 'Name your model here — the overlay is click-through, so you can type.'},
            {target: runButton, position: 'top', extra: statusBar,
              content: 'An element target works too, and `extra` punches a second hole over the status text.'},
            {target: 'noSuchControl', content: 'You never see this step: its target does not exist, so the ' +
              'tour skips forward and warns once in the console.'},
            {target: 'runButton', advanceOn: hasRun,
              content: 'Click Run — the tour advances by itself when the signal turns true. NEXT still works.'},
          ],
          onFinish: (r) => result.value = r,
        });
      }, {primary: true}),
      button('Finish now', () => tour?.finish()),
    ]),
    readout('onFinish', result, 'tour-result'),
    readout('hasRun', computed(() => String(hasRun.value))),
  ]);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'A running tour owns one scope for itself and one per step: the keyboard listener, the ' +
    'floating-ui autoUpdate and the advanceOn effect all die with them, so finishing the tour ' +
    'returns live scopes to the page baseline.'));
  main.append(button('Dispose sections', () => {
    tour?.finish();
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
