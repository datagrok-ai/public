/* The writable half of the panel (D-9): a form per section over the node's property model, every
   field committing through the engine, and the rows refreshed in place on every patch — a rebuilt
   row would drop the input the user is typing in (D-4). */
import * as grok from 'datagrok-api/grok';
import {Scope} from '../../core/scope.js';
import {div, h3, span} from '../../core/elements.js';
import {Accordion} from '../../components/containers/accordion.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import {PropertyGrid} from '../../components/forms/property-grid.js';
import type {PropDescriptor} from '../../components/forms/property-grid.js';
import type {PropertyLike} from '../../core/property-like.js';
import {propertyForm} from '../forms/object-form.js';
import type {FieldOverride, ObjectForm} from '../forms/object-form.js';
import type {SpecEditor, SpecPatch} from '../../spec/editor.js';
import type {SpecNodeRef} from './node-ref.js';
import {refreshPanel} from './handler.js';
import {lookGrid} from './look-grid.js';
import type {LookGrid} from './look-grid.js';
import {bindPicker, bindPickerButton} from './bind-picker.js';
import {eventEntry, eventPick, funcPicker, paramProps} from './func-picker.js';
import {bindRowsOf, bindsOf, commitOnChange, eventsOf, missingTable, paramBinds, paramValuesOf,
  paramsOf, propertyTier, propsFor, shownCommand, stringProps, unboundOf} from './prop-model.js';

/** {@link PropertyGrid.same}'s comparison — element-wise where both values are lists, identity
 * otherwise — for values the document may carry as arrays. */
const LISTY: PropDescriptor = {name: '', type: 'string_list'};
/** The automation handle of the "Add binding…" row and its picker button. */
const ADD_BINDING = 'add-binding';
/** Under the grid in Run mode, where an edit changes the live viewer and nothing else (R-c). */
const RUN_HINT = 'Edits here change the running view only and are not saved to the form — ' +
  'switch to Design to edit the form.';
/** What a bindable prop's name reads as in the "Add binding…" list: `xColumnName` → `X Column Name`. */
const words = (name: string): string => name.replace(/([a-z])([A-Z])/g, '$1 $2').replace(/^./, (c) => c.toUpperCase());

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

/** What {@link propEditors} answers: the sections, and the Bindings section drawn again from the
 * document — its rows ARE the document's bindings, so a pick adds one. */
export interface PanelEditors {
  sections: HTMLElement[];
  refresh: () => void;
}

/** The writable panel: a form per section, each field committing through the engine — and a
 * refusal putting the field back and saying why, as every canvas gesture does. */
export function propEditors(x: SpecNodeRef, editor: SpecEditor, events: Record<string, string>,
  panel: Scope): PanelEditors {
  const node = x.node;
  const owner = Scope.ambient!;
  const channels: FormChannel[] = [];
  const sections: HTMLElement[] = [];
  // the wiring sections fold: one row per bindable prop is a lot of panel for a node that binds
  // nothing, and it is what pushed the events off the bottom of the platform's pane
  let folds: Accordion | undefined;
  const channelOf = (props: PropertyLike[], values: Record<string, unknown>,
    read: () => Record<string, unknown>, patch: FormChannel['patch'],
    overrides?: Record<string, FieldOverride>): FormChannel => {
    // built before the form: a commit fired during construction must find its channel
    const channel = {values, read, patch} as FormChannel;
    channel.form = propertyForm(props, values, {condensed: true, overrides,
      onChanged: (name, value) => commit(editor, channel, name, value)});
    // the deterministic automation handle: every field is addressable by its prop name
    for (const prop of props)
      channel.form.input(prop.name)?.root.setAttribute('data-u2-prop', prop.name);
    channels.push(channel);
    return channel;
  };
  const place = (title: string, content: HTMLElement, expanded?: boolean): void => {
    if (expanded === undefined) {
      sections.push(h3(title), content);
      return;
    }
    if (folds === undefined) {
      folds = new Accordion();
      sections.push(folds.root);
    }
    // eagerly, never through the lazy builder: the fields refresh in place on every patch, so
    // they must exist whether or not the pane was ever opened
    folds.addPane(title, content, expanded);
  };
  const add = (title: string, props: PropertyLike[], values: Record<string, unknown>,
    read: () => Record<string, unknown>, patch: FormChannel['patch'],
    overrides?: Record<string, FieldOverride>, expanded?: boolean): FormChannel => {
    const channel = channelOf(props, values, read, patch, overrides);
    place(title, channel.form.root, expanded);
    return channel;
  };

  // a property-tier node's look is edited in the platform's own grid (VP-21): one grid commit is
  // one undo entry, its members pre-checked so a refused one is dropped alone, with a warning;
  // the commit runs from a microtask, so what the apply throws is owned here
  const tier = propertyTier(x);
  let look: LookGrid | undefined;
  if (tier) {
    look = lookGrid(x, panel, (patches) => {
      try {
        const accepted: SpecPatch[] = [];
        for (const patch of patches) {
          const refusal = editor.canApply(patch);
          if (refusal === null)
            accepted.push(patch);
          else
            grok.shell.warning(refusal);
        }
        editor.applyAll(accepted);
      } catch (e) {
        grok.shell.error(e instanceof Error ? e.message : String(e));
      }
    });
    if (x.instance.designTime)
      sections.push(h3('Properties'), look.root);
    else
      sections.push(h3('Properties (live)'), span(RUN_HINT, 'u2-designer-hint'), look.root);
  } else {
    for (const section of propsFor(x)) {
      add(section.title, section.props, section.values, section.read, (name, value) => {
        const current = node.props?.[name];
        // a component that reports '' for a prop the spec never carried would write that noise
        // into the document the moment the empty field is touched
        if (PropertyGrid.same(LISTY, current, value) || (value === '' && current === undefined))
          return null;
        return {op: 'set-prop', node, name, value};
      });
    }
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
      bindPickerButton(input, key, () => Scope.runWith(owner, () => bindPicker(x.instance,
        (path) => reissue(editor, {op: 'set-bind', node, name: key, path}, x))));
    }
  }

  // Bindings and events are guarded fields: a partial value ('c', 'cm', …) is a refusal, so
  // they commit on blur/Enter — per keystroke, the guard would fire on every prefix typed
  const bindings = (): FormChannel => {
    const bind = bindRowsOf(x);
    const rows = Object.keys(bind);
    // the rows are the section's for its lifetime: a cleared one stays, empty, until the next render
    const read = (): Record<string, unknown> => {
      const all = bindsOf(x);
      const shown: Record<string, unknown> = {};
      for (const name of rows)
        shown[name] = all[name] ?? '';
      return shown;
    };
    const channel = channelOf(stringProps(bind, 'Bound context path'), bind, read, (name, value) => {
      const path = value as string;
      return path === (node.bind?.[name] ?? '') ? null :
        {op: 'set-bind', node, name, path: path === '' ? undefined : path};
    }, commitOnChange(bind));
    // the path field is the power path; the picker beside it is the one that needs no grammar
    for (const name of rows) {
      const input = channel.form.input(name);
      if (input !== undefined) {
        bindPickerButton(input, name, () => Scope.runWith(owner, () =>
          bindPicker(x.instance, (path) => commit(editor, channel, name, path))));
      }
    }
    if (!tier && rows.length === 0)
      channel.form.addElement(span('Bind a parameter with the … beside it', 'u2-designer-hint'));
    // a property-tier node lists what it binds, not its forty bindable props: the next one is
    // picked here, and the row it lands on becomes the binding's — the section is drawn again
    if (tier) {
      const pick = new ChoiceInput({label: 'Add binding…', nullable: true,
        items: unboundOf(x).map((name) => ({value: name, label: words(name)})),
        tooltipText: 'The property to bind; … picks what it follows'});
      pick.root.setAttribute('data-u2-prop', ADD_BINDING);
      channel.form.add(pick);
      bindPickerButton(pick, ADD_BINDING, () => {
        const name = pick.value.peek();
        if (name === null) {
          grok.shell.warning('Add binding: pick the property first');
          return;
        }
        Scope.runWith(owner, () => bindPicker(x.instance,
          (path) => reissue(editor, {op: 'set-bind', node, name, path}, x)));
      });
    }
    return channel;
  };
  // the section is the document's bindings: built again in place when a pick adds one — under
  // its own scope, so the form it replaces is released with it
  let bindScope: Scope | undefined;
  let bound: FormChannel | undefined;
  panel.own(() => bindScope?.dispose());
  const host = div([], 'u2-designer-bindings');
  const refresh = (): void => {
    if (bound !== undefined)
      channels.splice(channels.indexOf(bound), 1);
    bindScope?.dispose();
    bindScope = new Scope();
    bound = Scope.runWith(bindScope, bindings);
    host.replaceChildren(bound.form.root);
  };
  // a source's params bind there too, so its section is up before the first one is picked
  if (Object.keys(bindRowsOf(x)).length > 0 || tier || params !== null) {
    refresh();
    // open where there is wiring to see — or where the one missing binding is what broke the node —
    // folded where every row is an empty invitation
    place('Bindings', host, Object.keys(node.bind ?? {}).length > 0 || missingTable(x));
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
    // before the node gate: a patch on the source re-renders the viewer too, and the grid would
    // otherwise edit a dead look
    look?.refresh();
    if (!applied || applied.structural || !applied.patches.some((p) => p.node === x.node))
      return;
    for (const channel of channels) {
      Object.assign(channel.values, channel.read());
      channel.form.refresh();
    }
  });
  return {sections, refresh: () => {
    refresh();
    look?.refresh();
  }};
}

function commit(editor: SpecEditor, channel: FormChannel, name: string, value: unknown): void {
  const patch = channel.patch(name, value);
  if (patch !== null)
    apply(editor, patch);
}

/** A patch the row it lands on becomes the binding's to show: the panel is drawn again in place —
 * the platform renders a `shell.o` write one gesture behind, so the write alone left the new row
 * unseen — and then asked for again, unless the selection already did that (a bind that flipped
 * the node between broken and built): one `shell.o` write per gesture. */
function reissue(editor: SpecEditor, patch: SpecPatch, x: SpecNodeRef): void {
  const shown = grok.shell.o;
  if (!apply(editor, patch))
    return;
  refreshPanel(x);
  if (grok.shell.o === shown)
    grok.shell.o = x;
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
