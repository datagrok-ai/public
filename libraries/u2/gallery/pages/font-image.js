import {signal, computed, Scope, Control} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {FontInput} from '../../src/components/inputs/font-input.js';
import {ImageInput} from '../../src/components/inputs/image-input.js';
import {TextInput} from '../../src/components/inputs/text-input.js';

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
injectOnce('u2-font-css', '../../css/font.css');

// a tiny inline SVG, so the page has a picture without reaching for the network
const SAMPLE = 'data:image/svg+xml,' + encodeURIComponent(
  '<svg xmlns="http://www.w3.org/2000/svg" width="120" height="80">' +
  '<rect width="120" height="80" fill="#2083D5"/>' +
  '<circle cx="60" cy="40" r="24" fill="#FFFFFF"/></svg>');

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
  main.append(el('h1', null, 'Font & image inputs'));
  const intro = el('p');
  intro.innerHTML = 'Two value editors ported from the platform. <code>FontInput</code> edits the ' +
    'strict font string viewer looks carry — <code>"&lt;weight&gt; &lt;style&gt; &lt;size&gt;px ' +
    '&lt;family&gt;"</code> — through a size box with presets, B / I toggles and a family select; ' +
    'generic families are written lowercase and unquoted. <code>ImageInput</code> previews the ' +
    'URL or data URL it holds, with a REMOVE link on hover. The password and auto-resize ' +
    '<code>TextInput</code> variants sit at the bottom.';
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

  const font = signal('bold normal 14px "Roboto"');
  section('Font', () => {
    const input = new FontInput({label: 'Font', bind: font});
    const sample = el('div', null, 'The quick brown fox jumps over the lazy dog');
    Scope.ambient.effect(() => sample.style.font = font.value);
    return [input, sample, readout('font', font)];
  });

  section('Font — custom families and sizes', () => {
    const input = new FontInput({label: 'Chart title',
      value: 'normal italic 18px "Open Sans"',
      families: ['Open Sans', 'Roboto Mono', 'Georgia'],
      sizes: [12, 18, 24, 36]});
    return [input, readout('value', computed(() => input.value.value))];
  });

  section('Font — a family the select does not offer', () => {
    const input = new FontInput({label: 'Legacy look', value: 'normal normal 11px "Comic Sans MS"'});
    return [input, readout('value', computed(() => input.value.value)),
      el('p', 'u2-gallery-status',
        'The unknown family joins the select as a hidden option, so it survives the round-trip.')];
  });

  const picture = signal(SAMPLE);
  section('Image', () => {
    const input = new ImageInput({label: 'Preview', bind: picture, alt: 'Sample'});
    return [input, readout('value length', computed(() => String(picture.value.length)))];
  });
  main.append(el('p', 'u2-gallery-status',
    'REMOVE appears on hover and clears the value; an empty value hides the preview.'));
  main.append(button('Restore the picture', () => picture.value = SAMPLE));

  section('Text — password', () => {
    const input = new TextInput({label: 'Password', password: true, value: 'hunter2'});
    return [input, readout('value', computed(() => input.value.value)),
      el('p', 'u2-gallery-status',
        'The eye toggles masking. The platform PASSWORD_PLACEHOLDER dance is a credentials ' +
        'convention and stays out of the core input.')];
  });

  section('Text — autoResize', () => {
    const input = new TextInput({label: 'Name', autoResize: true, placeholder: 'Compound',
      value: 'Aspirin'});
    return [input, el('p', 'u2-gallery-status',
      'The field grows with the text between minWidth (100) and maxWidth (300).')];
  });

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Each section is a Control.build(...) owner: disposing it releases the listeners, the size ' +
    'popup and the value effects, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
