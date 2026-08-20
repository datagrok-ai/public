/* The whole property-panel mechanism of the designer (DD1): selecting a node makes it the shell's
   current object, and this handler renders it — so spec nodes reach the context panel the way every
   other platform object does, with no bespoke panel plumbing. Read-only at P2. */
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {divV, h3, span} from '../../core/elements.js';
import {tableFromMap} from '../../components/table.js';
import type {PropertyLike} from '../../core/property-like.js';
import {SpecNodeRef} from './node-ref.js';

/** Heading for the properties a component declares without a category of their own. */
const UNGROUPED = 'Properties';

export class SpecNodeHandler extends DG.ObjectHandler<SpecNodeRef> {
  get type(): string {
    return 'u2-spec-node';
  }

  isApplicable(x: any): boolean {
    return x instanceof SpecNodeRef;
  }

  getCaption(x: SpecNodeRef): string {
    const caption = x.node.name ? `${x.node.name} (${x.node.tag})` : x.node.tag;
    return x.error() === null ? caption : `${caption} — broken`;
  }

  renderMarkup(x: SpecNodeRef): HTMLElement {
    return span(this.getCaption(x), 'u2-designer-caption');
  }

  renderProperties(x: SpecNodeRef): HTMLElement {
    const identity: Record<string, unknown> = {Tag: x.node.tag, Name: x.node.name ?? '', Path: x.path()};
    const error = x.error();
    if (error !== null)
      identity.Error = error;

    // the panel carries its own header: the platform shows one only for entities it knows
    const sections: HTMLElement[] = [this.renderMarkup(x), h3('Node'), tableFromMap(identity)];
    for (const [category, props] of SpecNodeHandler._grouped(x)) {
      const values: Record<string, unknown> = {};
      for (const prop of props)
        values[prop.caption ?? prop.friendlyName ?? prop.name] = SpecNodeHandler._text(prop.get?.(null));
      sections.push(h3(category), tableFromMap(values));
    }

    const events = SpecNodeHandler._events(x);
    if (Object.keys(events).length > 0)
      sections.push(h3('Events'), tableFromMap(events));
    const bind = x.node.bind;
    if (bind && Object.keys(bind).length > 0)
      sections.push(h3('Bindings'), tableFromMap({...bind}));
    return divV(sections, 'u2-designer-properties');
  }

  /** Declaration order, categories in the order they first appear. */
  private static _grouped(x: SpecNodeRef): Map<string, PropertyLike[]> {
    const built = x.built();
    const groups = new Map<string, PropertyLike[]>();
    if (!(built instanceof Control))
      return groups;
    for (const prop of built.getProperties()) {
      const key = prop.category ?? UNGROUPED;
      const group = groups.get(key);
      if (group)
        group.push(prop);
      else
        groups.set(key, [prop]);
    }
    return groups;
  }

  /** Every event the component declares, wired or not, plus anything the node wires beyond them. */
  private static _events(x: SpecNodeRef): Record<string, string> {
    const on = x.node.on ?? {};
    const events: Record<string, string> = {};
    for (const name of x.meta()?.events ?? [])
      events[name] = on[name] ?? '';
    for (const name of Object.keys(on))
      events[name] = on[name];
    return events;
  }

  private static _text(value: unknown): string {
    if (value === null || value === undefined)
      return '';
    return typeof value === 'object' ? JSON.stringify(value) : String(value);
  }
}

let registered = false;

/** Once per session: the platform keeps every handler it is given, and reopening the designer
 * would otherwise stack duplicates. */
export function registerSpecNodeHandler(): void {
  if (registered)
    return;
  registered = true;
  DG.ObjectHandler.register(new SpecNodeHandler());
}
