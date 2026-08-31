/* The writable half of the panel (D-9): a form per section over the node's property model, every
   field committing through the engine, and the rows refreshed in place on every patch — a rebuilt
   row would drop the input the user is typing in (D-4). */
import * as grok from 'datagrok-api/grok';
import {Scope} from '../../core/scope.js';
import {div, h3, span} from '../../core/elements.js';
import {Accordion} from '../../components/containers/accordion.js';
import type {AccordionPane} from '../../components/containers/accordion.js';
import {ChoiceInput} from '../../components/inputs/choice-input.js';
import {PropertyGrid} from '../../components/forms/property-grid.js';
import type {PropDescriptor} from '../../components/forms/property-grid.js';
import type {NamedProperty} from '../../core/widget-like.js';
import {propertyForm} from '../forms/object-form.js';
import type {FieldOverride, ObjectForm} from '../forms/object-form.js';
import {APPEARANCE_CATEGORY} from '../../spec/appearance.js';
import type {SpecEditor, SpecPatch} from '../../spec/editor.js';
import type {SpecNodeRef} from './node-ref.js';
import {refreshPanel} from './handler.js';
import {lookGrid} from './look-grid.js';
import type {LookGrid} from './look-grid.js';
import {bindPicker, bindPickerButton} from './bind-picker.js';
import {BindChip} from './bind-chip.js';
import {eventEntry, eventPick, funcPicker, paramProps} from './func-picker.js';
import {bindsOf, commitOnChange, eventsOf, missingTable, paramBinds, paramValuesOf,
  paramsOf, propertyTier, propsFor, sharedAppearance, shownCommand, stringProps, unboundOf} from './prop-model.js';

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
  let badge: (() => void) | undefined;
  const channelOf = (props: NamedProperty[], values: Record<string, unknown>,
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
  const place = (title: string, content: HTMLElement, expanded?: boolean): AccordionPane | undefined => {
    if (expanded === undefined) {
      sections.push(h3(title), content);
      return undefined;
    }
    if (folds === undefined) {
      folds = new Accordion();
      sections.push(folds.root);
    }
    // eagerly, never through the lazy builder: the fields refresh in place on every patch, so
    // they must exist whether or not the pane was ever opened
    return folds.addPane(title, content, expanded);
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
  }

  // every property row is visually bindable (UB-7a): the picker button rides the rail, and a
  // bound row renders the chip in place of its editor (UB-7b). Bind/unbind flips a row between
  // chip and editor, and the platform ignores a same-identity `shell.o` write — so the rows are
  // rebuilt here, in place, under their own scope, whenever a set-bind lands on the node
  const pickFor = (name: string, label = name) => () => Scope.runWith(owner, () => bindPicker(x.instance,
    (path) => reissue(editor, {op: 'set-bind', node, name, path}, x),
    {prop: label, current: node.bind?.[name]}));
  const chipsFor = (props: NamedProperty[]): Record<string, FieldOverride> | undefined => {
    let overrides: Record<string, FieldOverride> | undefined;
    for (const prop of props) {
      const name = prop.name;
      const path = node.bind?.[name];
      if (path === undefined)
        continue;
      const chip = new BindChip({path, label: prop.caption ?? prop.friendlyName ?? name,
        twoWay: (prop as {twoWay?: boolean}).twoWay === true,
        onPick: pickFor(name),
        onClear: () => reissue(editor, {op: 'set-bind', node, name, path: undefined}, x)});
      chip.clearButton.dataset.u2Unbind = name;
      (overrides ??= {})[name] = {input: chip};
    }
    return overrides;
  };
  const rowButtons = (channel: FormChannel, props: NamedProperty[]): void => {
    for (const prop of props) {
      // a subBindable head binds through its dotted parameter rows only — a button here would
      // dead-end at the canApply refusal
      const meta = prop as {subBindable?: boolean, bindable?: boolean};
      if (meta.subBindable === true && meta.bindable !== true)
        continue;
      const input = channel.form.input(prop.name);
      if (input !== undefined)
        bindPickerButton(input, prop.name, pickFor(prop.name));
    }
  };
  // shared-prop identity, not the section title: a component-own prop that declares the shared
  // category keeps ordinary ''-writing semantics
  const shared = sharedAppearance(x);
  let rowScope: Scope | undefined;
  let rowChannels: FormChannel[] = [];
  panel.own(() => rowScope?.dispose());
  // the h3/pane placement happens once — the titles are stable across rebuilds (propsFor groups
  // by category + Parent tag); every rebuild swaps the form inside its stable host
  const hosts = new Map<string, HTMLElement>();
  if (!tier) {
    for (const section of propsFor(x, true)) {
      const h = div([], 'u2-designer-rows');
      hosts.set(section.title, h);
      if (section.title !== APPEARANCE_CATEGORY) {
        place(section.title, h);
        continue;
      }
      // the shared group folds, collapsed until something is assigned — or a component-own prop
      // joined the section by declaring its category; the badge counts the assigned shared values
      const props = section.props;
      const count = (): number => props.filter((p) => shared.has(p.name) &&
        (node.props?.[p.name] !== undefined || node.bind?.[p.name] !== undefined)).length;
      const hasOwn = props.some((p) => !shared.has(p.name));
      const pane = place(section.title, h, count() > 0 || hasOwn);
      if (pane !== undefined) {
        const label = pane.root.querySelector('.u2-accordion-title') as HTMLElement;
        badge = () => {
          const n = count();
          label.textContent = n > 0 ? `${APPEARANCE_CATEGORY} (${n})` : APPEARANCE_CATEGORY;
        };
        badge();
      }
    }
  }
  let paramsHost: HTMLElement | undefined;
  if (paramsOf(x) !== null)
    place('Parameters', paramsHost = div([], 'u2-designer-rows'));
  const buildRows = (): void => {
    for (const channel of rowChannels)
      channels.splice(channels.indexOf(channel), 1);
    rowChannels = [];
    rowScope?.dispose();
    rowScope = new Scope();
    Scope.runWith(rowScope, () => {
      if (!tier) {
        for (const section of propsFor(x, true)) {
          const patch: FormChannel['patch'] = (name, value) => {
            // a bound prop is the context's to change: a stale row never writes a literal over it
            if (node.bind?.[name] !== undefined)
              return null;
            const current = node.props?.[name];
            // clearing a shared appearance field DELETES the prop (Ruling 8): absent means platform styling
            if (shared.has(name) && (value === '' || value === null))
              return current === undefined ? null : {op: 'set-prop', node, name, value: undefined};
            // a component that reports '' for a prop the spec never carried would write that noise
            // into the document the moment the empty field is touched
            if (PropertyGrid.same(LISTY, current, value) || (value === '' && current === undefined))
              return null;
            return {op: 'set-prop', node, name, value};
          };
          const channel = channelOf(section.props, section.values, section.read, patch,
            chipsFor(section.props));
          rowChannels.push(channel);
          rowButtons(channel, section.props);
          hosts.get(section.title)!.replaceChildren(channel.form.root);
        }
      }
      // what the function this source names is called with: one typed editor per declared param,
      // and the same picker the props use for the ones that follow a value instead
      const params = paramsOf(x);
      if (params !== null) {
        const values = paramValuesOf(node, params.prop);
        const channel = channelOf(
          paramProps(params.inputs, values, paramBinds(node, params.prop)), values,
          () => paramValuesOf(node, params.prop), (name, value) => {
            const current = paramValuesOf(node, params.prop);
            if (PropertyGrid.same(LISTY, current[name], value) ||
              (value === '' && current[name] === undefined))
              return null;
            return {op: 'set-prop', node, name: params.prop, value: {...current, [name]: value}};
          });
        rowChannels.push(channel);
        for (const prop of params.inputs) {
          const input = channel.form.input(prop.name);
          const key = `${params.prop}.${prop.name}`;
          if (input !== undefined)
            bindPickerButton(input, key, pickFor(key, prop.name));
        }
        paramsHost?.replaceChildren(channel.form.root);
      }
    });
  };
  buildRows();

  // Bindings and events are guarded fields: a partial value ('c', 'cm', …) is a refusal, so
  // they commit on blur/Enter — per keystroke, the guard would fire on every prefix typed.
  // One presentation for every node (UB-7c): bound rows only — free-text, the power path for
  // re-pointing — plus one "Add binding…" over every unbound declared prop
  const bindings = (): FormChannel => {
    const bind = bindsOf(x);
    const rows = Object.keys(bind);
    const offered = unboundOf(x);
    // the rows are the section's for its lifetime: a cleared one stays, empty, until the next render
    const read = (): Record<string, unknown> => {
      const all = bindsOf(x);
      const shown: Record<string, unknown> = {};
      for (const name of rows)
        shown[name] = all[name] ?? '';
      return shown;
    };
    // the same humanized labels the "Add binding…" list speaks; a dotted sub-bind key keeps its
    // parameter name visible past the head. The prop key stays the row's identity untouched
    const caption = (name: string): string => {
      const dot = name.indexOf('.');
      return dot < 0 ? words(name) : `${words(name.slice(0, dot))} › ${name.slice(dot + 1)}`;
    };
    const channel = channelOf(
      stringProps(bind, 'Bound context path').map((p) => ({...p, caption: caption(p.name)})),
      bind, read, (name, value) => {
        const path = value as string;
        return path === (node.bind?.[name] ?? '') ? null :
          {op: 'set-bind', node, name, path: path === '' ? undefined : path};
      }, commitOnChange(bind));
    // the path field is the power path; the picker beside it is the one that needs no grammar
    for (const name of rows) {
      const input = channel.form.input(name);
      if (input !== undefined) {
        bindPickerButton(input, name, () => Scope.runWith(owner, () =>
          bindPicker(x.instance, (path) => commit(editor, channel, name, path),
            {prop: name, current: node.bind?.[name]})));
      }
    }
    // the next binding is picked here, and the row it lands on becomes the binding's — the
    // section is drawn again
    if (tier || offered.length > 0) {
      const pick = new ChoiceInput({label: 'Add binding…', nullable: true,
        items: offered.map((name) => ({value: name, label: words(name)})),
        tooltipText: 'The property to bind; … picks what it follows'});
      pick.root.setAttribute('data-u2-prop', ADD_BINDING);
      channel.form.add(pick);
      bindPickerButton(pick, ADD_BINDING, () => {
        const name = pick.value.peek();
        if (name === null) {
          grok.shell.warning('Add binding: pick the property first');
          return;
        }
        pickFor(name)();
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
  const hasBinds = Object.keys(bindsOf(x)).length > 0;
  if (hasBinds || tier || paramsOf(x) !== null || unboundOf(x).length > 0) {
    refresh();
    // open where there is wiring to see — or where the one missing binding is what broke the node —
    // folded where every row is an empty invitation
    place('Bindings', host, hasBinds || missingTable(x));
  }

  if (Object.keys(events).length > 0) {
    const channel = channelOf(
      stringProps(events,
        'The function this event calls — pick one with the … button, or type cmd: plus its name'),
      events,
      () => eventsOf(x), (name, value) => {
        const command = value as string;
        return command === shownCommand(node.on?.[name]) ? null :
          {op: 'set-event', node, event: name, command: command === '' ? undefined : command};
      }, commitOnChange(events));
    // always open: this is the section the folded one above was hiding
    place('Events', channel.form.root, true);
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
    // a bind or unbind flips rows between chip and editor: rebuild them — never mid-typing, since
    // typing commits set-prop/set-event, which refresh in place below (D-4)
    if (applied.patches.some((p) => p.node === x.node && p.op === 'set-bind'))
      buildRows();
    for (const channel of channels) {
      Object.assign(channel.values, channel.read());
      channel.form.refresh();
    }
    // the count badge follows set/clear from any surface; the fold itself is never forced
    badge?.();
  });
  return {sections, refresh: () => {
    buildRows();
    refresh();
    look?.refresh();
  }};
}

/** False on a refusal — what keeps the picker that produced the value open for a correction. */
function commit(editor: SpecEditor, channel: FormChannel, name: string, value: unknown): boolean {
  const patch = channel.patch(name, value);
  return patch === null || apply(editor, patch);
}

/** A patch that changes which rows are chips: the document changes, then the panel redraws itself
 * in place (`refreshPanel` → the row rebuild and the Bindings section). The platform ignores a
 * same-identity `shell.o` write, so the panel never relies on one. False on a refusal. */
function reissue(editor: SpecEditor, patch: SpecPatch, x: SpecNodeRef): boolean {
  if (!apply(editor, patch))
    return false;
  refreshPanel(x);
  return true;
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
