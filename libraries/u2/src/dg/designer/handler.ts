/* The whole property-panel mechanism of the designer (DD1): selecting a node makes it the shell's
   current object, and this handler renders it — so spec nodes reach the context panel the way every
   other platform object does, with no bespoke panel plumbing. The panel writes exactly when the ref
   carries an editor (D-9): the absence of that channel is the read-only mode, and both paths read
   the same property model, so a plain HTML node shows its props either way. */
import * as DG from 'datagrok-api/dg';
import {Component} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {divV, h3, span} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {tableFromMap} from '../../components/collections/table.js';
import {SpecNodeRef, SpecNodesRef} from './node-ref.js';
import {propEditors} from './prop-editors.js';
import {eventsOf, propsFor} from './prop-model.js';
import type {PropSection} from './prop-model.js';
import {sourceStatus} from './source-status.js';

export {propsFor} from './prop-model.js';
export type {PropSection} from './prop-model.js';

/** The live panel's editors and their effects: one per session, replaced on every render. The
 * platform gives no teardown hook, so the designer view owns the last one's disposal. */
let panel: Scope | undefined;

export function disposePanel(): void {
  panel?.dispose();
  panel = undefined;
}

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
    disposePanel();
    const identity: Record<string, unknown> = {Tag: x.node.tag, Name: x.node.name ?? '', Path: x.path()};
    const error = x.error();
    if (error !== null)
      identity.Error = error;

    // the panel carries its own header: the platform shows one only for entities it knows
    const sections: HTMLElement[] = [this.renderMarkup(x), h3('Node'), tableFromMap(identity)];
    const model = propsFor(x);
    const events = eventsOf(x);
    const editor = x.editor;
    if (editor === undefined) {
      sections.push(...SpecNodeHandler._status(x, undefined));
      sections.push(...SpecNodeHandler._tables(model, events, {...x.node.bind}));
    } else {
      const scope = new Scope();
      panel = scope;
      sections.push(...SpecNodeHandler._status(x, scope));
      sections.push(...Scope.runWith(scope, () => propEditors(x, editor, model, events, scope)));
    }
    return divV(sections, 'u2-designer-properties');
  }

  /** A data source's pulse: the one LIVE section of the panel, because a source runs on its own
   * clock and the canvas shows the same empty fields whether it is loading, was never run or
   * failed. Read-only and cheap — a redrawn table touches no editor, so the in-place refresh the
   * form fields rely on is untouched. */
  private static _status(x: SpecNodeRef, scope: Scope | undefined): HTMLElement[] {
    const built = x.built();
    if (!(built instanceof Component) || sourceStatus(built) === null)
      return [];
    // a source with no `designData` prop never runs anything the user did not ask for: saying so
    // is what keeps its ABSENCE from reading as a missing feature
    const policy = x.meta()?.props.some((prop) => prop.name === 'designData') === true ? null :
      'always live — this source has no design-time policy';
    const host = divV([], 'u2-designer-status');
    const draw = (): void => {
      const at = sourceStatus(built)!;
      const shown: Record<string, unknown> = {};
      if (at.state !== '')
        shown.State = at.state === 'idle' ? 'idle — not run yet; use Refresh' : at.state;
      if (at.count !== null)
        shown[at.counts === 'rows' ? 'Rows' : 'Entities'] = at.count;
      if (at.error !== null)
        shown.Error = at.error;
      if (policy !== null)
        shown['Design data'] = policy;
      host.textContent = '';
      host.append(tableFromMap(shown));
    };
    if (scope === undefined)
      draw();
    else
      scope.effect(draw);
    return [h3('Status'), host];
  }

  /** The read-only panel, built from the same model — which is what finally gives a plain HTML node
   * (and a node that failed to build) properties to show. */
  private static _tables(model: PropSection[], events: Record<string, string>,
    bind: Record<string, string>): HTMLElement[] {
    const sections: HTMLElement[] = [];
    for (const section of model) {
      const values: Record<string, unknown> = {};
      for (const prop of section.props)
        values[prop.name] = text(section.values[prop.name]);
      sections.push(h3(section.title), tableFromMap(values));
    }
    if (Object.keys(events).length > 0)
      sections.push(h3('Events'), tableFromMap(events));
    if (Object.keys(bind).length > 0)
      sections.push(h3('Bindings'), tableFromMap(bind));
    return sections;
  }
}

/** The multiple-objects convention over a multi-selection: how many, of what, and where —
 * read-only, no property editing across a multi-selection. */
export class SpecNodesHandler extends DG.ObjectHandler<SpecNodesRef> {
  get type(): string {
    return 'u2-spec-nodes';
  }

  isApplicable(x: any): boolean {
    return x instanceof SpecNodesRef;
  }

  getCaption(x: SpecNodesRef): string {
    return `${x.refs.length} nodes`;
  }

  renderMarkup(x: SpecNodesRef): HTMLElement {
    return span(this.getCaption(x), 'u2-designer-caption');
  }

  renderProperties(x: SpecNodesRef): HTMLElement {
    // the multi panel replaces the single-node one — its editors must not outlive it
    disposePanel();
    const counts: Record<string, unknown> = {};
    for (const ref of x.refs)
      counts[ref.node.tag] = Number(counts[ref.node.tag] ?? 0) + 1;
    return divV([
      this.renderMarkup(x),
      h3('Tags'), tableFromMap(counts),
      h3('Nodes'), divV(x.refs.map((ref) => span(ref.path()))),
    ], 'u2-designer-properties');
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
  DG.ObjectHandler.register(new SpecNodesHandler());
}
