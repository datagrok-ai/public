/* The toolbox palette: every component of the instance's registry as an accordion — one pane per
   `category`, filtered by one search box across panes. A row carries its tag in `data-u2-tag` and
   nothing else — the drag layer lives in the view, so the palette stays a plain, testable list
   (a mousedown on a pane header finds no `.u2-palette-item` and falls through to the toggle). */
import {Control} from '../../core/component.js';
import {div, span} from '../../core/elements.js';
import {Accordion} from '../../components/containers/accordion.js';
import type {AccordionPane} from '../../components/containers/accordion.js';
import {TextInput} from '../../components/inputs/text-input.js';
import type {ComponentMeta, Registry} from '../../spec/registry.js';

const UNGROUPED = 'Other';
/** The panes in the order the toolbox shows them: what lays out, what shows data, what edits it;
 * a category not named here follows, first-seen. */
const CATEGORY_ORDER = ['Containers', 'Viewers', 'Inputs', 'Actions', 'Display', 'Data'];
/** The one pane sorted by label: it lists every platform viewer, in whatever order the platform
 * registered them. */
const SORTED = 'Viewers';

/** Collapse memory, session-level by design (F9): survives a palette rebuild within the session,
 * gone with the page — the durable medium waits for P5's asset. */
const expandedByCategory = new Map<string, boolean>();

interface Group {
  category: string;
  pane: AccordionPane;
  items: {el: HTMLElement, text: string}[];
}

export class Palette extends Control {
  readonly filter: TextInput;

  private readonly _groups: Group[] = [];

  constructor(reg: Registry) {
    super();
    this.root.classList.add('u2-palette');
    this.root.dataset.u2 = 'palette';
    this.filter = this.runInScope(() => new TextInput({search: true, inline: true, placeholder: 'Filter'}));
    const accordion = this.runInScope(() => new Accordion());
    accordion.root.classList.add('u2-palette-list');
    this.root.append(this.filter.root, accordion.root);

    for (const [category, metas] of Palette._grouped(reg)) {
      const items = metas.map((meta) => ({el: Palette._item(meta),
        text: `${meta.tag} ${Palette._label(meta)} ${meta.description} ${meta.usage ?? ''}`.toLowerCase()}));
      const pane = accordion.addPane(category, div(items.map((item) => item.el)),
        expandedByCategory.get(category) ?? true);
      this._groups.push({category, pane, items});
      // a toggle is remembered only while no filter holds panes open — a force-expanded pane
      // is the filter's doing, not the user's choice
      this.effect(() => {
        const on = pane.expanded.value;
        if (this.filter.value.peek().trim() === '')
          expandedByCategory.set(category, on);
      });
    }
    this.effect(() => this._apply(this.filter.value.value.trim().toLowerCase()));
  }

  private _apply(query: string): void {
    const filtering = query !== '';
    for (const group of this._groups) {
      let shown = 0;
      for (const item of group.items) {
        const hidden = filtering && !item.text.includes(query);
        item.el.hidden = hidden;
        if (!hidden)
          shown++;
      }
      // a filter whose matches sit inside a collapsed pane would be useless: matching panes
      // force-open while the query holds, and clearing it restores the remembered state
      group.pane.root.hidden = filtering && shown === 0;
      if (filtering && shown > 0)
        group.pane.expanded.value = true;
      else if (!filtering)
        group.pane.expanded.value = expandedByCategory.get(group.category) ?? true;
    }
  }

  /** Registration order within a category, the categories in {@link CATEGORY_ORDER} and the
   * unknown ones after, in the order they first appear. */
  private static _grouped(reg: Registry): Map<string, ComponentMeta[]> {
    const groups = new Map<string, ComponentMeta[]>();
    for (const meta of reg.metas()) {
      const key = meta.category ?? UNGROUPED;
      const group = groups.get(key);
      if (group)
        group.push(meta);
      else
        groups.set(key, [meta]);
    }
    const rank = (category: string): number => {
      const at = CATEGORY_ORDER.indexOf(category);
      return at === -1 ? CATEGORY_ORDER.length : at;
    };
    const ordered = new Map([...groups].sort((a, b) => rank(a[0]) - rank(b[0])));
    ordered.get(SORTED)?.sort((a, b) => Palette._label(a).localeCompare(Palette._label(b)));
    return ordered;
  }

  /** What the row says: the meta's own label where it carries one (a viewer's platform name),
   * the tag suffix otherwise. */
  private static _label(meta: ComponentMeta): string {
    return meta.label ?? (meta.tag.startsWith('u2-') ? meta.tag.slice(3) : meta.tag);
  }

  private static _item(meta: ComponentMeta): HTMLElement {
    const el = div([], 'u2-palette-item');
    const glyph = meta.icon?.();
    if (glyph)
      el.append(div([glyph], 'u2-palette-item-icon'));
    el.append(span(Palette._label(meta), 'u2-palette-item-label'));
    el.title = meta.usage ?? meta.description;
    el.dataset.u2Tag = meta.tag;
    return el;
  }
}
