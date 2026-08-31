import * as grok from 'datagrok-api/grok';
import {SpecContext, signal} from '@datagrok-libraries/u2';
import {platformContext} from '@datagrok-libraries/u2/src/dg/index.js';

/* What the designer opens on: one representative dg-ui/1 form — a splitter over a named form, a
   prop bound two-way to context data, a command-wired button and a text area — small enough to read
   whole in the "Open…" dialog, and varied enough that the structure tree, the property panel and the
   hit-test all have something to show. Every interactive node carries a `name`: that is what the
   inspector selects by. */
export const DESIGNER_SPEC = {
  $schema: 'dg-ui/1',
  root: {
    tag: 'u2-splitter',
    name: 'layout',
    props: {direction: 'horizontal', sizes: [0.45, 0.55]},
    children: [
      {tag: 'u2-panel', name: 'left', children: [
        {tag: 'h2', props: {text: 'Compound'}},
        {tag: 'u2-form', name: 'details', children: [
          {tag: 'u2-text-input', name: 'nameInput', props: {label: 'Name'},
            bind: {value: '$.reagent'}},
          {tag: 'u2-number-input', name: 'doseInput',
            props: {label: 'Dose', value: 250, min: 0, max: 1000, postfix: 'mg'}},
          {tag: 'u2-choice-input', name: 'stageInput',
            props: {label: 'Stage', items: ['Discovery', 'Phase I', 'Phase II'], value: 'Phase I'}},
          {tag: 'u2-bool-input', name: 'activeInput', props: {label: 'Active', value: true, switch: true}},
        ]},
        {tag: 'u2-button', name: 'saveButton', props: {text: 'Save', primary: true},
          on: {click: 'cmd:save'}},
      ]},
      {tag: 'u2-panel', name: 'right', children: [
        {tag: 'h2', props: {text: 'Notes'}},
        {tag: 'u2-text-area', name: 'notesInput',
          props: {label: 'Notes', placeholder: 'Dissolved in saline, stored at 4 °C.', autoGrow: true}},
        {tag: 'u2-text-area', name: 'runLog',
          props: {label: 'Run log', placeholder: 'Commands log their runs here.', autoGrow: true},
          bind: {value: '$.demoRunLog'}},
      ]},
    ],
  },
};

/* The balloon autohides in seconds; this is the durable trace of every command and every function
   the designer fires. It outlives one designer session on purpose — `u2Record` (package.ts) writes
   here from a spec the user wired themselves, and the app is closed and reopened while doing it. */
const demoRunLog = signal('');

export function appendRunLog(line: string): void {
  const shown = `${new Date().toLocaleTimeString()} ${line}`;
  demoRunLog.value = `${demoRunLog.peek()}${demoRunLog.peek() === '' ? '' : '\n'}${shown}`;
}

/** What the spec above resolves against: the `$.` data the name field binds two-way to, and the
 * `cmd:` names it wires — those run in Run mode, and are inert under the glass pane. Platform-backed,
 * so a button the function picker wires reaches `grok.functions.call`. A brand-new form still offers
 * these keys — they belong to the app, not to the form — so the picker shows them under App data and
 * the run log carries `demo` in its name rather than reading as the designer's own plumbing. */
export function designerContext(): SpecContext {
  const reagent = signal('Aspirin');
  return platformContext({
    data: {reagent, demoRunLog},
    commands: {save: () => {
      const message = `cmd:save ran for "${reagent.peek()}" — the button carries the name, not the code.`;
      grok.shell.info(message);
      appendRunLog(message);
    }},
  });
}
