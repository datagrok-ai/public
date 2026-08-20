import * as grok from 'datagrok-api/grok';
import {SpecContext, signal} from '@datagrok-libraries/u2';

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
      ]},
    ],
  },
};

/** What the spec above resolves against: the `$.` data the name field binds two-way to, and the
 * `cmd:` names it wires — those run in Run mode, and are inert under the glass pane. */
export function designerContext(): SpecContext {
  const reagent = signal('Aspirin');
  return new SpecContext({
    data: {reagent},
    commands: {save: () => grok.shell.info(`cmd:save ran for "${reagent.peek()}" — ` +
      'the button carries the name, not the code.')},
  });
}
