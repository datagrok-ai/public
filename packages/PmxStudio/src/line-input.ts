/* The tumor model, picked from the atlas's models: a field showing the chosen one that opens, in place, the models
   grouped by class (Syngeneic, CDX, PDX), filtered by class and found by name, tumor type or organ, each with the
   papers and animals its growth prior comes from. The value is the line's name, as the functions' `line` param. */
import * as grok from 'datagrok-api/grok';
import {ButtonGroup, div, Input, InputOptions, Reason, signal, span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {TumorModel, tumorModels} from './priors';

export type ModelClass = TumorModel['class'];
export const CLASSES: {cls: ModelClass, label: string, tip: string}[] = [
  {cls: 'syngeneic', label: 'Syngeneic', tip: 'A mouse tumor line in immunocompetent mice'},
  {cls: 'cdx', label: 'CDX', tip: 'A human cell line in immunodeficient mice'},
  {cls: 'pdx', label: 'PDX', tip: 'A patient\'s tumor in immunodeficient mice'},
];

/** The models of class `cls` (null: all) with every word of `query` in their name, tumor type or organ, by class in
 *  CLASSES' order, each class in the atlas's order. */
export function modelsFound(models: {[line: string]: TumorModel}, cls: ModelClass | null, query: string):
  {cls: ModelClass, lines: string[]}[] {
  const words = query.toLowerCase().split(/\s+/).filter((w) => w);
  const found = (line: string, m: TumorModel) => {
    const text = `${line} ${m.tumour_type} ${m.organ}`.toLowerCase();
    return words.every((w) => text.includes(w));
  };
  return CLASSES.filter((c) => cls === null || c.cls === cls)
    .map((c) => ({cls: c.cls, lines: Object.entries(models).filter(([l, m]) => m.class === c.cls && found(l, m))
      .map(([l]) => l)}))
    .filter((g) => g.lines.length > 0);
}

const LAST = 'lastLine';
/** The model the user picked last (Plan opens on it), or null. */
export const lastLine = (): string | null => grok.userSettings.getValue('PmxStudio', LAST) || null;

/** The tumors of `line` as a protocol names them: syngeneic tumors, or cell-line or patient-derived xenografts. */
export const tumorsOf = (line: string, m?: TumorModel) => !m ? `${line} tumors` : m.class === 'syngeneic' ?
  `${line} syngeneic tumors` : `${line} ${m.class === 'cdx' ? 'cell-line' : 'patient-derived'} xenografts ` +
  `(${m.class.toUpperCase()}, ${m.tumour_type})`;

export class LineInput extends Input<string | null, InputOptions<string | null>> {
  private declare _models: ReturnType<typeof signal<{[line: string]: TumorModel} | null>>;
  /** The atlas's models, once read. */
  get models(): ReadonlySignal<{[line: string]: TumorModel} | null | undefined> {return this._models;}

  constructor(options: InputOptions<string | null>) {
    super(options, 'MC38');
    this.root.dataset.u2 = 'pmx-line-input';
    void tumorModels().then((m) => this._models.value = m, () => this._models.value = {});
  }

  protected createEditor(): HTMLElement {
    // (createEditor runs inside the base constructor, before this class's fields are set)
    this._models = signal<{[line: string]: TumorModel} | null>(null);
    const open = signal(false);
    const cls = signal<ModelClass | null>(null);
    const query = signal('');
    // the line the keyboard is on in the open list
    const active = signal<string | null>(null);

    const field = document.createElement('button');
    field.type = 'button';
    field.className = 'pmx-line-field';
    field.setAttribute('aria-haspopup', 'listbox');
    const name = span('', 'pmx-line-name');
    const kind = span('', 'pmx-line-kind');
    field.append(name, kind, span('▾', 'pmx-line-chevron'));

    const groups = new ButtonGroup({items: [{id: 'all', text: 'All'}, ...CLASSES.map((c) =>
      ({id: c.cls, text: c.label, tooltip: c.tip}))], toggle: 'single'});
    groups.selected.value = ['all'];
    const search = document.createElement('input');
    search.type = 'search';
    search.className = 'pmx-line-search';
    search.placeholder = 'Search: name, tumor type, organ';
    search.setAttribute('aria-label', 'Search tumor models');
    const list = div([], 'pmx-line-list');
    list.setAttribute('role', 'listbox');
    const panel = div([groups.root, search, list], 'pmx-line-panel');

    const visible = () => modelsFound(this._models.peek() ?? {}, cls.peek(), query.peek()).flatMap((g) => g.lines);
    const pick = (line: string) => {
      this.value = line;
      grok.userSettings.add('PmxStudio', LAST, line);
      void grok.userSettings.flush();
      open.value = false;
      field.focus();
    };
    const move = (by: number) => {
      const lines = visible();
      if (lines.length === 0)
        return;
      const at = lines.indexOf(active.peek() ?? '');
      active.value = lines[Math.min(lines.length - 1, Math.max(0, at < 0 ? 0 : at + by))];
      list.querySelector(`[data-line="${CSS.escape(active.peek()!)}"]`)?.scrollIntoView({block: 'nearest'});
    };

    this.runInScope(() => {
      this.scope.listen(field, 'click', () => open.value = !open.peek());
      this.scope.listen(search, 'input', () => query.value = search.value);
      this.scope.listen(panel, 'keydown', (ev) => {
        const e = ev as KeyboardEvent;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          move(e.key === 'ArrowDown' ? 1 : -1);
        } else if (e.key === 'Enter' && active.peek()) {
          e.preventDefault();
          e.stopPropagation();
          pick(active.peek()!);
        } else if (e.key === 'Escape') {
          e.stopPropagation();
          open.value = false;
          field.focus();
        }
      });
      this.effect(() => {
        const g = groups.selected.value[0];
        cls.value = !g || g === 'all' ? null : g as ModelClass;
      });
      this.effect(() => {
        const v = this.value;
        const m = this._models.value?.[v ?? ''];
        name.textContent = v ?? '';
        kind.textContent = m ? `${m.tumour_type} · ${CLASSES.find((c) => c.cls === m.class)?.label ?? m.class}` : '';
      });
      this.effect(() => {
        panel.hidden = !open.value;
        field.setAttribute('aria-expanded', String(open.value));
        if (open.value) {
          search.value = '';
          query.value = '';
          active.value = this.peek();
          search.focus();
        }
      });
      // the class buttons with their counts; one with no models stays, off, saying so
      this.effect(() => {
        const models = this._models.value ?? {};
        for (const c of CLASSES) {
          const b = groups.root.querySelector<HTMLElement>(`[data-id="${c.cls}"]`);
          const n = Object.values(models).filter((m) => m.class === c.cls).length;
          const label = b?.querySelector('.u2-btn-group-label');
          if (!b || !label)
            continue;
          label.textContent = `${c.label} ${n}`;
          Reason.apply(b, this._models.value && n === 0 ? `No ${c.label} models in the atlas yet` : null);
        }
      });
      this.effect(() => {
        const models = this._models.value;
        const at = active.value;
        const chosen = this.value;
        if (!models) {
          list.replaceChildren(span('Loading the atlas’s models…', 'pmx-line-empty'));
          return;
        }
        const found = modelsFound(models, cls.value, query.value);
        list.replaceChildren(...found.length === 0 ? [span('No model matches', 'pmx-line-empty')] :
          found.flatMap((g) => [span(CLASSES.find((c) => c.cls === g.cls)!.label, 'pmx-line-group'),
            ...g.lines.map((line) => {
              const m = models[line];
              const item = div([div([span(line, 'pmx-line-item-name'), span(m.tumour_type, 'pmx-line-item-type')]),
                span(`${m.papers} papers · ${m.animals} animals`, 'pmx-line-item-evidence')], 'pmx-line-item');
              item.dataset.line = line;
              item.setAttribute('role', 'option');
              item.setAttribute('aria-selected', String(line === chosen));
              item.classList.toggle('is-active', line === at);
              item.title = `${line}: ${m.tumour_type} (${m.organ}), ${m.host_class} host. Its growth prior comes ` +
                `from ${m.papers} papers, ${m.animals} animals.`;
              item.addEventListener('click', () => pick(line));
              return item;
            })]));
      });
    });
    return div([field, panel], 'pmx-line-input');
  }
}
