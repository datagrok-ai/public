import {signal, computed, Scope, Component} from '../../src/index.js';
import {divV, divH, span, button} from '../../src/core/elements.js';
import {TextInput} from '../../src/components/text-input.js';
import {Form} from '../../src/components/form.js';
import {Wizard} from '../../src/components/wizard.js';

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
injectOnce('u2-form-css', '../../css/form.css');
injectOnce('u2-wizard-css', '../../css/wizard.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function paragraph(text) {
  const p = el('div');
  p.style.padding = '0 0 var(--dg-space-m)';
  p.textContent = text;
  return p;
}

/** The same four steps both inline and in a dialog: free text, a gated form, lazy content,
 * and a summary that reads the inputs' signals. */
function makeWizard(logLine, outcome) {
  const notes = new TextInput({label: 'Notes', placeholder: 'Anything at all'});
  const email = new TextInput({label: 'Email', placeholder: 'name@example.com'});
  email.addValidator((v) => /^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(v) ? null : 'Enter a valid email address');
  const contact = new Form().addAll([email]);
  const flavor = signal('(not built yet)');

  const steps = [
    {
      id: 'notes',
      title: 'Notes',
      content: divV([paragraph('Nothing gates this step — NEXT is always enabled.'), notes.root]),
    },
    {
      id: 'contact',
      title: 'Contact',
      content: divV([paragraph('NEXT is gated by the form validity signal.'), contact.root]),
      canProceed: contact.validity,
    },
    {
      id: 'options',
      title: 'Options',
      content: () => {
        flavor.value = `built at ${new Date().toLocaleTimeString()}`;
        logLine(`options step built (${flavor.value}) — never rebuilt afterwards`);
        return divV([paragraph('Built lazily on first visit; go back and forth, the log stays at one line.'),
          span(computed(() => `content ${flavor.value}`))]);
      },
      onActivate: (step) => logLine(`onActivate: ${step.id}`),
    },
    {
      id: 'summary',
      title: 'Summary',
      content: divV([
        paragraph('Reads the signals of the inputs from the earlier steps.'),
        span(computed(() => `Notes: ${notes.value.value || '(empty)'}`)),
        span(computed(() => `Email: ${email.value.value || '(empty)'}`)),
        span(computed(() => `Options: ${flavor.value}`)),
      ]),
    },
  ];

  return new Wizard({
    steps,
    onFinish: () => outcome.value = `finished with ${email.value.peek() || '(no email)'}`,
    onCancel: () => outcome.value = 'cancelled',
  });
}

export async function render(main) {
  main.append(el('h1', null, 'Wizard'));
  const intro = el('p');
  intro.innerHTML = 'Linear steps with a gate per step: <code>canProceed</code> takes a signal or a ' +
    'callback and is normalized to one blocking reason, which disables NEXT and shows itself next to ' +
    'the buttons. Panels are hidden rather than detached, so what you type survives navigation, and ' +
    'lazy steps are built once. <code>Enter</code> inside the content advances (<code>Ctrl+Enter</code> ' +
    'finishes on the last step); ←/→ move between the step markers, and a completed marker is clickable.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const log = el('div', 'u2-gallery-status');
  const logLine = (text) => log.append(el('div', null, text));

  const parts = [];
  const section = (title, builder) => {
    main.append(el('h2', null, title));
    const component = Component.build(builder);
    parts.push(component);
    main.append(component.root);
    return component;
  };

  section('Inline', () => {
    const outcome = signal('(nothing yet)');
    const wizard = makeWizard(logLine, outcome);
    wizard.root.style.height = '320px';
    wizard.root.style.maxWidth = '620px';
    wizard.root.style.border = 'var(--dg-border-width) solid var(--dg-border)';
    return [
      wizard,
      readout('currentStep', wizard.currentStep),
      readout('completed', computed(() => String(wizard.completed.value))),
      readout('outcome', outcome),
    ];
  });

  section('In a dialog', () => {
    const outcome = signal('(nothing yet)');
    const wizard = makeWizard(logLine, outcome);
    return [
      divH([button('Open the same wizard in a dialog',
        () => wizard.openInDialog('New compound', {width: 520}))]),
      el('p', 'u2-gallery-status',
        'The dialog shows no buttons of its own: BACK / NEXT / FINISH and CANCEL are the wizard\'s ' +
        'footer. Esc, the ✕ and CANCEL all report a cancel; FINISH closes the dialog. Reopening ' +
        'returns to the step you left.'),
      readout('currentStep', wizard.currentStep),
      readout('completed', computed(() => String(wizard.completed.value))),
      readout('outcome', outcome),
    ];
  });

  main.append(el('h2', null, 'Lazy build log'));
  main.append(log);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Both wizards were built inside a Component.build(...) builder, so the sections own them: ' +
    'disposing releases the per-step scopes (including whatever a lazy step built), the rail and ' +
    'content listeners and the dialog — live scopes return to the page baseline.'));
  main.append(button('Dispose sections', () => {
    for (const part of parts)
      part.dispose();
    refresh();
  }));
  refresh();
}
