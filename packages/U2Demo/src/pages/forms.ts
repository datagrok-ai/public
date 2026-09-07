import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  signal, computed,
  divV, span, h3, button, bigButton,
  TextInput, BoolInput, NumberInput, ChoiceInput,
  Form, PropertyGrid,
} from '@datagrok-libraries/u2';
import {propertyForm, IProperty} from '@datagrok-libraries/u2/src/dg/index.js';
import {readout} from './common';

export function formPage(): HTMLElement {
  const first = new TextInput({label: 'First name', name: 'first', value: 'Ada'});
  first.addValidator((v) => v.trim() ? null : 'Required');
  const last = new TextInput({label: 'Last name', name: 'last'});
  last.addValidator((v) => v.trim() ? null : 'Required');
  const age = new NumberInput({label: 'Age', name: 'age', mode: 'int', min: 0, max: 120, value: 36});
  const role = new ChoiceInput({label: 'Role', name: 'role', items: ['Viewer', 'Editor', 'Admin'], value: 'Editor'});
  const subscribed = new BoolInput({label: 'Subscribe', name: 'subscribed', switch: true, value: true});

  const form = new Form();
  form.addAll([first, last, age, role, subscribed]);
  form.addButtons((row) => {
    row.append(bigButton('Submit', () => {
      if (form.validate())
        grok.shell.info(JSON.stringify(form.getValues()));
    }));
    row.append(button('Reset', () => form.setValues({first: 'Ada', last: '', age: 36, role: 'Editor',
      subscribed: true})));
  });

  return divV([
    h3('Form'),
    span('Aligned label column, aggregated validity, Enter submits.', 'u2demo-hint'),
    form,
    readout('form validity', computed(() => form.validity.value ?? 'valid')),
  ], 'u2demo-page');
}

export function propertyGridPage(): HTMLElement {
  const pg = new PropertyGrid();
  pg.setProperties([
    {name: 'title', type: 'string', description: 'Viewer title'},
    {name: 'opacity', type: 'double', min: 0, max: 1},
    {name: 'width', type: 'int', min: 100, max: 2000, category: 'Layout'},
    {name: 'height', type: 'int', min: 100, max: 2000, category: 'Layout'},
    {name: 'showLegend', type: 'bool', category: 'Legend'},
    {name: 'position', type: 'choice', choices: ['left', 'right', 'top', 'bottom'], category: 'Legend'},
  ], {title: 'Scatter plot', opacity: 0.8, width: 600, height: 400, showLegend: true, position: 'right'});

  return divV([
    h3('PropertyGrid'),
    span('The value record is replaced on every edit — bind to it directly.', 'u2demo-hint'),
    pg,
    readout('last change', computed(() => {
      const c = pg.onChanged.value;
      return c ? `${c.name} → ${String(c.value)}` : '(none)';
    })),
    readout('values', computed(() => JSON.stringify(pg.values.value))),
  ], 'u2demo-page');
}

export function objectFormPage(): HTMLElement {
  const group = DG.Group.create(`u2demo-${Math.random().toString(36).slice(2, 8)}`);
  const props: IProperty[] = [
    {name: 'friendlyName', caption: 'Name', type: 'string', nullable: false,
      get: (g) => (g as DG.Group).friendlyName, set: (g, v) => (g as DG.Group).friendlyName = v as string},
    {name: 'personal', caption: 'Personal', type: 'bool',
      get: (g) => (g as DG.Group).personal, set: (g, v) => (g as DG.Group).personal = v as boolean},
    {name: 'hidden', caption: 'Hidden', type: 'bool',
      get: (g) => (g as DG.Group).hidden, set: (g, v) => (g as DG.Group).hidden = v as boolean},
  ];
  const form = propertyForm(props, group);
  const result = signal('(not saved)');
  form.addButtons((row) => {
    row.append(
      bigButton('SAVE', async () => {
        if (!form.validate()) {
          result.value = 'Fix validation first';
          return;
        }
        result.value = 'Saving…';
        try {
          const saved = await grok.dapi.groups.save(group);
          result.value = `Saved: ${saved.friendlyName} (${saved.id})`;
        } catch (e) {
          result.value = `Save failed: ${e instanceof Error ? e.message : e}`;
        }
      }),
      button('Delete', async () => {
        try {
          await grok.dapi.groups.delete(group);
          result.value = 'Deleted';
        } catch (e) {
          result.value = `Delete failed: ${e instanceof Error ? e.message : e}`;
        }
      }));
  });
  return divV([
    span('propertyForm(): everything below is generated from IProperty metadata over a real DG.Group — ' +
      'editors chosen by type, required from nullable: false, edits write through to the entity, ' +
      'SAVE persists it via grok.dapi.groups.', 'u2demo-hint'),
    span('SAVE and Delete need a server: in local mode dapi has nobody to talk to, so both fail — ' +
      'the reason lands in the result line below. Open the app on a stand to persist for real.',
    'u2demo-hint'),
    form.root,
    readout('result', result),
  ], 'u2demo-page');
}
