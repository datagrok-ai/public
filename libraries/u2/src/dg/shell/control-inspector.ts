/* Click a u2 control, see its properties: the DOM → component door (`Control.forElement`) meets
   the platform's object channel — an ObjectHandler renders a two-way PropertyGrid in the context
   panel. Opt-in: a plugin that wants the inspector calls registerControlInspector() once. */
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {divV, h3, span} from '../../core/elements.js';
import {PropertyGrid, PropDescriptor} from '../../components/forms/property-grid.js';
import {registry} from '../../spec/registry.js';
import {registerAll} from '../../spec/registrations.js';

const EDITABLE = new Set(['string', 'int', 'double', 'bool', 'string_list']);

/** The PropertyGrid rows a control answers with — the gallery inspector's mapping made reusable.
 * A control built by hand rather than through `registry.create` still gets the registry's
 * metadata, but a prop backed by a private field or by a constructor option nobody passed reads
 * nothing: a row that can neither show a value nor take one is worse than no row at all. */
export function controlPropDescriptors(control: Control): PropDescriptor[] {
  control.componentMeta ??= registry.get(`u2-${control.root.dataset.u2 ?? ''}`);
  const target = control.propertyTarget;
  return control.getProperties()
    .filter((p) => (EDITABLE.has(p.type ?? '') || p.choices) && p.get?.(target) !== undefined)
    .map((p) => ({name: p.name, type: (p.choices ? 'choice' : p.type) as PropDescriptor['type'],
      choices: p.choices, min: p.min, max: p.max, category: p.category,
      description: p.description, readonly: p.set == null}));
}

/** A grid over one control, two-way where the property has a setter. Disposed by the caller. */
export function controlProperties(control: Control): PropertyGrid {
  const descriptors = controlPropDescriptors(control);
  const props = control.getProperties();
  const grid = new PropertyGrid({onChanged: (change) =>
    props.find((p) => p.name === change.name)?.set?.(control.propertyTarget, change.value)});
  grid.setProperties(descriptors, Object.fromEntries(props.map((p) => [p.name, p.get?.(control.propertyTarget)])));
  return grid;
}

/** The control under `target`, bounded by `within` — `Control.forElement` walks to the document
 * root, so an unbounded click on the whitespace between two controls resolves to whatever shell
 * control happens to be an ancestor (an app's own splitter, its view host). Nothing outside the
 * subtree, and never the boundary itself. */
export function controlAt(target: Element | null, within: Element): Control | undefined {
  // npm datagrok-api typings lag the repo (forElement ships in 1.28); the runtime is DG.U2's
  const control = (Control as any).forElement(target) as Control | undefined;
  return control && control.root !== within && within.contains(control.root) ? control : undefined;
}

class NoControl {}

/** What a caller makes current when it wants the panel to say the object under the pointer is not
 * a u2 component. An inspector click that resolves to nothing is better off doing nothing —
 * replacing whatever the panel holds with this is destructive, not informative. */
export const noControl = new NoControl();

let panel: Scope | undefined;
/** Live {@link registerControlInspector} registrations; the platform cannot drop a handler, so
 * with none left the handlers stop claiming anything. */
let users = 0;

export function disposePanel(): void {
  panel?.dispose();
  panel = undefined;
}

class ControlHandler extends DG.ObjectHandler<Control> {
  get type(): string { return 'u2-control'; }

  // u2 controls only: every DG.Widget is a Control too, and a viewer or a platform widget keeps
  // the platform's own handler. Disjoint from the designer's SpecNodeRef and U2Demo's PropRow
  isApplicable(x: any): boolean { return users > 0 && Control.is(x) && !(x instanceof DG.Widget); }

  getCaption(x: Control): string { return x.componentMeta?.tag ?? x.root.dataset.u2 ?? 'control'; }

  renderProperties(x: Control): HTMLElement {
    disposePanel();
    // the header goes AFTER this call: reading the properties is what stamps `componentMeta` from
    // the registry, and the caption reads it — built first, the first render loses the `u2-` tag
    const descriptors = controlPropDescriptors(x);
    const caption = this.getCaption(x);
    // a control the registry does not describe (VirtualList, VirtualTree, Combobox — not spec
    // components) has no property metadata to read: say that, rather than render an empty grid
    if (descriptors.length === 0) {
      return divV([h3(caption), span(`No inspectable properties: ${caption} is not a registered ` +
        'u2 component, so nothing describes its properties.', 'u2-async-empty')], 'u2-inspector');
    }
    const scope = panel = new Scope(undefined, {label: 'inspector'});
    const grid = scope.build(() => controlProperties(x));
    const root = divV([h3(caption), grid], 'u2-inspector');
    // the context panel kills what it shows before it renders the next object
    DG.Widget.registerCleanup(root, () => scope.dispose());
    return root;
  }
}

class NoControlHandler extends DG.ObjectHandler<NoControl> {
  get type(): string { return 'u2-no-control'; }

  isApplicable(x: any): boolean { return users > 0 && x instanceof NoControl; }

  getCaption(): string { return 'Not a u2 component'; }

  renderProperties(): HTMLElement {
    disposePanel();
    return divV([span('This element is not a u2 component, so there is nothing to inspect. ' +
      'Click a u2 control to see its live properties here.', 'u2-async-empty')], 'u2-inspector');
  }
}

let registered = false;

/** Turns the inspector on; the handlers are registered with the platform once per session.
 * Answers the unregister: once every registration is undone, u2 controls go back to whatever
 * else handles them and the open panel is disposed. */
export function registerControlInspector(): () => void {
  if (!registered) {
    registered = true;
    registerAll();
    DG.ObjectHandler.register(new ControlHandler());
    DG.ObjectHandler.register(new NoControlHandler());
  }
  users++;
  let done = false;
  return () => {
    if (done)
      return;
    done = true;
    if (--users === 0)
      disposePanel();
  };
}
