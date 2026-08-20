/* The toolbox palette: every component of the instance's registry as an accordion — one pane per
   `category`, filtered by one search box across panes. A row carries its tag in `data-u2-tag` and
   nothing else — the drag layer lives in the view, so the palette stays a plain, testable list
   (a mousedown on a pane header finds no `.u2-palette-item` and falls through to the toggle). */
import {Control} from '../../core/component.js';
import {div} from '../../core/elements.js';
import {Accordion} from '../../components/containers/accordion.js';
import type {AccordionPane} from '../../components/containers/accordion.js';
import {TextInput} from '../../components/inputs/text-input.js';
import type {ComponentMeta, Registry} from '../../spec/registry.js';

const UNGROUPED = 'Other';

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
    this.filter = this.run(() => new TextInput({search: true, inline: true, placeholder: 'Filter'}));
    const accordion = this.run(() => new Accordion());
    accordion.root.classList.add('u2-palette-list');
    this.root.append(this.filter.root, accordion.root);

    for (const [category, metas] of Palette._grouped(reg)) {
      const items = metas.map((meta) => ({el: Palette._item(meta),
        text: `${meta.tag} ${meta.description}`.toLowerCase()}));
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

  /** Registration order within a category, categories in the order they first appear. */
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
    return groups;
  }

  private static _item(meta: ComponentMeta): HTMLElement {
    const el = div([], 'u2-palette-item');
    el.textContent = meta.tag.startsWith('u2-') ? meta.tag.slice(3) : meta.tag;
    el.title = meta.description;
    el.dataset.u2Tag = meta.tag;
    return el;
  }
}
