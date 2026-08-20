/* Handler-backed objects on u2 surfaces: `handlerRenderer` adapts the platform's
   ObjectHandler to the core `ObjectRenderer` seam, `chip` is the u2 counterpart of
   `ui.render(x)`, and `entityInput` composes both with `dapiSource` into the generic
   server-entity picker (OBJECTS.md in core/docs/features/ui2). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Control} from '../core/component.js';
import {Tooltip} from '../core/tooltip.js';
import {ObjectRenderer} from '../core/object-renderer.js';
import {TypeAhead} from '../components/inputs/typeahead.js';
import {dapiSource, DapiSourceLike, DapiSourceOptions} from './dapi-source.js';

/** An `ObjectRenderer` over the platform's ObjectHandler registry. The handler is resolved
 * once, from the first item seen — never `isApplicable` per item — so one renderer instance
 * serves one homogeneous collection. Handler output is foreign DOM: rendered by plugins,
 * never disposed or introspected here, and a throwing handler degrades to caption text. */
export class HandlerRenderer<T> implements ObjectRenderer<T> {
  private _handler: DG.ObjectHandler | null | undefined;

  handlerFor(x: T): DG.ObjectHandler | null {
    if (this._handler === undefined)
      this._handler = DG.ObjectHandler.forEntity(x);
    return this._handler;
  }

  caption(x: T): string {
    try {
      const handler = this.handlerFor(x);
      return handler ? handler.getCaption(x) : String(x);
    } catch (_) {
      return String(x);
    }
  }

  icon(x: T): HTMLElement { return this._render(x, (h) => h.renderIcon(x)); }
  listItem(x: T): HTMLElement { return this._render(x, (h) => h.renderListItem(x)); }
  markup(x: T): HTMLElement { return this._render(x, (h) => h.renderMarkup(x)); }
  tooltip(x: T): HTMLElement { return this._render(x, (h) => h.renderTooltip(x)); }
  card(x: T): HTMLElement { return this._render(x, (h) => h.renderCard(x)); }

  private _render(x: T, fn: (handler: DG.ObjectHandler) => HTMLElement): HTMLElement {
    try {
      const handler = this.handlerFor(x);
      if (handler)
        return fn(handler);
    } catch (_) { /* foreign plugin code; fall back to caption */ }
    const el = document.createElement('span');
    el.textContent = this.caption(x);
    return el;
  }
}

export function handlerRenderer<T>(): HandlerRenderer<T> {
  return new HandlerRenderer<T>();
}

export interface ChipOptions {
  renderer?: ObjectRenderer<any>;
  /** Click makes the object the shell's current object (context panel). Default true. */
  currentOnClick?: boolean;
}

/** An inline representation of a platform object — handler markup, handler tooltip through
 * the u2 tooltip service, click sets `grok.shell.o`, and a semantic automation id
 * (`data-entity`, plus `data-entity-id` for a `DG.Entity`). */
export class EntityChip extends Control {
  constructor(x: any, options: ChipOptions = {}) {
    super(document.createElement('span'));
    const renderer = options.renderer ?? handlerRenderer();
    this.root.classList.add('u2-chip');
    this.root.dataset.u2 = 'chip';
    this.root.dataset.entity = renderer instanceof HandlerRenderer ?
      renderer.handlerFor(x)?.type ?? typeof x : typeof x;
    if (x instanceof DG.Entity && x.id)
      this.root.dataset.entityId = x.id;
    this.root.append(renderer.markup ? renderer.markup(x) : renderer.caption(x));
    if (renderer.tooltip)
      Tooltip.bind(this.root, () => renderer.tooltip!(x), this.scope);
    if (options.currentOnClick ?? true) {
      const onClick = () => grok.shell.o = x;
      this.root.addEventListener('click', onClick);
      this.own(() => this.root.removeEventListener('click', onClick));
    }
  }
}

export function chip(x: any, options: ChipOptions = {}): EntityChip {
  return new EntityChip(x, options);
}

/** A card-sized preview of a platform object — the handler's `renderCard` (markup, then
 * caption, as fallbacks), with the same current-object click and automation ids as
 * {@link EntityChip}. The platform styles the card content itself. */
export class EntityCard extends Control {
  constructor(x: any, options: ChipOptions = {}) {
    super();
    const renderer = options.renderer ?? handlerRenderer();
    this.root.classList.add('u2-entity-card');
    this.root.dataset.u2 = 'entity-card';
    this.root.dataset.entity = renderer instanceof HandlerRenderer ?
      renderer.handlerFor(x)?.type ?? typeof x : typeof x;
    if (x instanceof DG.Entity && x.id)
      this.root.dataset.entityId = x.id;
    this.root.append(renderer.card ? renderer.card(x) :
      renderer.markup ? renderer.markup(x) : renderer.caption(x));
    if (options.currentOnClick ?? true) {
      const onClick = () => grok.shell.o = x;
      this.root.addEventListener('click', onClick);
      this.own(() => this.root.removeEventListener('click', onClick));
    }
  }
}

export function entityCard(x: any, options: ChipOptions = {}): EntityCard {
  return new EntityCard(x, options);
}

export interface EntityInputOptions extends DapiSourceOptions {
  minChars?: number;
  debounceMs?: number;
}

/** Picks a server entity by type-ahead: `entityInput('Group…', () => grok.dapi.groups)`.
 * Fetches through {@link dapiSource}, renders through the object's handler; `selected`
 * holds the entity itself. */
export function entityInput<T>(placeholder: string, factory: () => DapiSourceLike<T>,
  options: EntityInputOptions = {}): TypeAhead<T> {
  return new TypeAhead<T>({
    source: dapiSource(factory, options),
    renderer: handlerRenderer<T>(),
    placeholder,
    minChars: options.minChars,
    debounceMs: options.debounceMs,
  });
}
