/* The writable half of the panel (D-9): a form per section over the node's property model, every
   field committing through the engine, and the rows refreshed in place on every patch — a rebuilt
   row would drop the input the user is typing in (D-4). */
import * as grok from 'datagrok-api/grok';
import {Scope} from '../../core/scope.js';
import {h3} from '../../core/elements.js';
import {Accordion} from '../../components/containers/accordion.js';
import {PropertyGrid} from '../../components/forms/property-grid.js';
import type {PropDescriptor} from '../../components/forms/property-grid.js';
import type {PropertyLike} from '../../core/property-like.js';
import {propertyForm} from '../object-form.js';
import type {FieldOverride, ObjectForm} from '../object-form.js';
import type {SpecEditor, SpecPatch} from '../../spec/editor.js';
import type {SpecNodeRef} from './node-ref.js';
import {bindPicker, bindPickerButton} from './bind-picker.js';
import {eventEntry, eventPick, funcPicker, paramProps} from './func-picker.js';
import {bindsOf, commitOnChange, eventsOf, paramBinds, paramValuesOf, paramsOf, shownCommand,
  stringProps} from './prop-model.js';
import type {PropSection} from './prop-model.js';

/** {@link PropertyGrid.same}'s comparison — element-wise where both values are lists, identity
 * otherwise — for values the document may carry as arrays. */
const LISTY: PropDescriptor = {name: '', type: 'string_list'};

/** One editable section of the panel — the form, the snapshot behind it, and the patch a field
 * edit becomes. */
interface FormChannel {
  form: ObjectForm;
  values: Record<string, unknown>;
  /** What the document says now: where a refused edit goes back to, and what an undo refreshes to. */
  read: () => Record<string, unknown>;
  /** Null for an edit that is not worth a patch (a no-op, or noise the node never carried). */
  patch: (name: string, value: unknown) => SpecPatch | null;
}

/** The writable panel: a form per section, each field committing through the engine — and a
 * refusal putting the field back and saying why, as every canvas gesture does. */
export function propEditors(x: SpecNodeRef, editor: SpecEditor, model: PropSection[],
  events: Record<string, string>, panel: Scope): HTMLElement[] {
  const node = x.node;
  const owner = Scope.ambient!;
  const channels: FormChannel[] = [];
  const sections: HTMLElement[] = [];
  // the wiring sections fold: one row per bindable prop is a lot of panel for a node that binds
  // nothing, and it is what pushed the events off the bottom of the platform's pane
  let folds: Accordion | undefined;
  const add = (title: string, props: PropertyLike[], values: Record<string, unknown>,
    read: () => Record<string, unknown>, patch: FormChannel['patch'],
    overrides?: Record<string, FieldOverride>, expanded?: boolean): FormChannel => {
    // built before the form: a commit fired during construction must find its channel
    const channel = {values, read, patch} as FormChannel;
    channel.form = propertyForm(props, values, {condensed: true, overrides,
      onChanged: (name, value) => commit(editor, channel, name, value)});
    // the deterministic automation handle: every field is addressable by its prop name
    for (const prop of props)
      channel.form.input(prop.name)?.root.setAttribute('data-u2-prop', prop.name);
    channels.push(channel);
    if (expanded === undefined)
      sections.push(h3(title), channel.form.root);
    else {
      if (folds === undefined) {
        folds = new Accordion();
        sections.push(folds.root);
      }
      // eagerly, never through the lazy builder: the fields refresh in place on every patch, so
      // they must exist whether or not the pane was ever opened
      folds.addPane(title, channel.form.root, expanded);
    }
    return channel;
  };

  for (const section of model) {
    add(section.title, section.props, section.values, section.read, (name, value) => {
      const current = node.props?.[name];
      // a component that reports '' for a prop the spec never carried would write that noise
      // into the document the moment the empty field is touched
      if (PropertyGrid.same(LISTY, current, value) || (value === '' && current === undefined))
        return null;
      return {op: 'set-prop', node, name, value};
    });
  }

  // what the function this source names is called with: one typed editor per declared param, and
  // the same picker the props use for the ones that follow a value instead
  const params = paramsOf(x);
  if (params !== null) {
    const values = paramValuesOf(node, params.prop);
    const channel = add('Parameters',
      paramProps(params.inputs, values, paramBinds(node, params.prop)), values,
      () => paramValuesOf(node, params.prop), (name, value) => {
        const current = paramValuesOf(node, params.prop);
        if (PropertyGrid.same(LISTY, current[name], value) ||
          (value === '' && current[name] === undefined))
          return null;
        return {op: 'set-prop', node, name: params.prop, value: {...current, [name]: value}};
      });
    for (const prop of params.inputs) {
      const input = channel.form.input(prop.name);
      if (input === undefined)
        continue;
      const key = `${params.prop}.${prop.name}`;
      bindPickerButton(input, key, () => Scope.runWith(owner, () => bindPicker(x.instance, (path) => {
        // the row becomes the binding's to show, so the panel has to be built again from it
        if (apply(editor, {op: 'set-bind', node, name: key, path}))
          grok.shell.o = x;
      })));
    }
  }

  // Bindings and events are guarded fields: a partial value ('c', 'cm', …) is a refusal, so
  // they commit on blur/Enter — per keystroke, the guard would fire on every prefix typed
  const bind = bindsOf(x);
  if (Object.keys(bind).length > 0) {
    const channel = add('Bindings', stringProps(bind, 'Bound context path'), bind,
      () => bindsOf(x), (name, value) => {
        const path = value as string;
        return path === (node.bind?.[name] ?? '') ? null :
          {op: 'set-bind', node, name, path: path === '' ? undefined : path};
      // open where there is wiring to see, folded where every row is an empty invitation
      }, commitOnChange(bind), Object.keys(node.bind ?? {}).length > 0);
    // the path field is the power path; the picker beside it is the one that needs no grammar
    for (const name of Object.keys(bind)) {
      const input = channel.form.input(name);
      if (input !== undefined) {
        bindPickerButton(input, name, () => Scope.runWith(owner, () =>
          bindPicker(x.instance, (path) => commit(editor, channel, name, path))));
      }
    }
  }

  if (Object.keys(events).length > 0) {
    const channel = add('Events',
      stringProps(events,
        'The function this event calls — pick one with the … button, or type cmd: plus its name'),
      events,
      () => eventsOf(x), (name, value) => {
        const command = value as string;
        return command === shownCommand(node.on?.[name]) ? null :
          {op: 'set-event', node, event: name, command: command === '' ? undefined : command};
      // always open: this is the section the folded one above was hiding
      }, commitOnChange(events), true);
    // the field names a command of any tier by hand; the picker beside it names a platform
    // function and what to call it with
    for (const name of Object.keys(events)) {
      const input = channel.form.input(name);
      if (input !== undefined) {
        bindPickerButton(input, name, () => Scope.runWith(owner, () => funcPicker({
          title: `${name}: pick a function`, instance: x.instance, pick: eventPick(node.on?.[name]),
          onPick: (picked) => apply(editor,
            {op: 'set-event', node, event: name, command: eventEntry(picked)}),
        })), 'Pick a function');
      }
    }
  }

  // an undo, a redo or an edit made elsewhere refreshes the fields in place: rebuilding them
  // would drop the input the user is typing in (D-4)
  panel.effect(() => {
    const applied = editor.onDidApply.value;
    if (!applied || applied.structural || !applied.patches.some((p) => p.node === x.node))
      return;
    for (const channel of channels) {
      Object.assign(channel.values, channel.read());
      channel.form.refresh();
    }
  });
  return sections;
}

function commit(editor: SpecEditor, channel: FormChannel, name: string, value: unknown): void {
  const patch = channel.patch(name, value);
  if (patch !== null)
    apply(editor, patch);
}

/** The one funnel every panel edit takes, whichever field or picker produced it. */
function apply(editor: SpecEditor, patch: SpecPatch): boolean {
  const refusal = editor.canApply(patch);
  if (refusal === null) {
    editor.apply(patch);
    // clearing wiring by emptying the field is designed (D-9), but it must never be silent —
    // an unwired button in Run mode gives no feedback at all
    if (patch.op === 'set-event' && patch.command === undefined)
      grok.shell.info(`${patch.event}: command removed`);
    else if (patch.op === 'set-bind' && patch.path === undefined)
      grok.shell.info(`${patch.name}: binding removed`);
    return true;
  }
  // the field keeps what was typed: a refusal is an invitation to correct it, and reverting
  // would wipe the text under the user (the document was never touched — nothing to revert)
  grok.shell.warning(refusal);
  return false;
}
