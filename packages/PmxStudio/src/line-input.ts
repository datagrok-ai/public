/* The tumor model, picked from the atlas's models: a field showing the chosen one that opens, in place, the models
   grouped by class (Syngeneic, CDX, PDX), filtered by class and found by name, tumor type or organ, each with the
   papers and animals its growth prior comes from. The value is the line's name, as the functions' `line` param;
   last, a model not in the atlas, planned from its own history (history-dialog.ts), whose value is its key. */
import {ButtonGroup, div, icon, Input, InputOptions, Reason, signal, span, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';

import {historyDialog} from './history-dialog';
import {HISTORY_PREFIX, TumorModel, enginePriors, historyKeys, historyModelName, tumorModels} from './priors';
import {setUserSetting, userSetting} from './ui';

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
export const lastLine = (): string | null => {
  const line = userSetting(LAST) || null;
  // a model planned from its history is not kept for the next plan: its history lives in its plan file
  return line?.startsWith(HISTORY_PREFIX) ? null : line;
};

/** The tumors of `line` as a protocol names them: syngeneic tumors, or cell-line or patient-derived xenografts. */
export function tumorsOf(line: string, m?: TumorModel): string {
  const own = historyModelName(line);
  if (own)
    return `${own} patient-derived xenografts (PDX, planned from its own history)`;
  if (!m)
    return `${line} tumors`;
  return m.class === 'syngeneic' ? `${line} syngeneic tumors` : `${line} ${m.class === 'cdx' ? 'cell-line' :
    'patient-derived'} xenografts (${m.class.toUpperCase()}, ${m.tumour_type})`;
}

/** The models whose ranges are a little less certain, said beside the model's name (R7 freeze record, limitation 3:
 *  the T1 fit-stability test not passed). */
export const FIT_NOTES: {[line: string]: string} = {
  '4T1': 'The data for this model fit almost equally well with slightly different growth curves, so its ranges ' +
    'are a little less certain than for other models.',
  'LLC': 'One part of this model\'s fit sits at its limit, so its ranges are a little less certain than for other ' +
    'models.',
};

/** Pickers made so far: each one's element ids its own. */
let pickers = 0;

export class LineInput extends Input<string | null, InputOptions<string | null>> {
  private declare _models: ReturnType<typeof signal<{[line: string]: TumorModel} | null>>;
  /** The atlas's models, once read. */
  get models(): ReadonlySignal<{[line: string]: TumorModel} | null | undefined> {return this._models;}

  constructor(options: InputOptions<string | null>) {
    // no model assumed: a study read from a file, or one without a plan, waits for one picked (QA 0.1.3 round 3)
    super(options, null);
    this.root.dataset.u2 = 'pmx-line-input';
    // each model with the papers and animals of its growth prior, the ones the plan's protocol cites (the models'
    // own listing counts a narrower set for some PDX models: 47 papers against the prior's 48)
    void Promise.all([tumorModels(), enginePriors()]).then(([m, p]) => this._models.value = Object.fromEntries(
      Object.entries(m).map(([l, x]) => [l, {...x, papers: p.lines[l]?.papers ?? x.papers,
        animals: p.lines[l]?.animals ?? x.animals}])), () => this._models.value = {});
    // a model the atlas does not have (a plan file edited by hand, a model since removed) is not kept as if it were
    // one: the field says so and the plan waits for a model from the list (QA 0.1.3 final designer #30)
    this.addValidator((v) => {
      const models = this._models.value;
      return !v || !models || !Object.keys(models).length || v in models || v.startsWith(HISTORY_PREFIX) ? null :
        `"${v}" is not a tumor model of the atlas: pick one from the list`;
    });
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
    // the same grey mark as Plan's lab mark (plan-view labMark), on the models FIT_NOTES names
    const note = icon('info-circle', {cls: 'pmx-line-note'});
    note.addEventListener('click', (e) => e.stopPropagation());
    field.append(name, note, kind, span('▾', 'pmx-line-chevron'));

    const groups = this.build(() => new ButtonGroup({items: [{id: 'all', text: 'All'}, ...CLASSES.map((c) =>
      ({id: c.cls, text: c.label, tooltip: c.tip}))], toggle: 'single'}));
    groups.selected.value = ['all'];
    const search = document.createElement('input');
    search.type = 'search';
    search.className = 'pmx-line-search';
    search.placeholder = 'Search: name, tumor type, organ';
    search.setAttribute('aria-label', 'Search tumor models');
    const list = div([], 'pmx-line-list');
    list.setAttribute('role', 'listbox');
    list.setAttribute('aria-label', 'Tumor models');
    // focus stays in the search while the arrows move: the option they are on is named to a screen reader from there
    // (ids of this picker's own: two pickers on a page do not share them)
    const ids = `pmx-line-${++pickers}`;
    list.id = `${ids}-list`;
    search.setAttribute('aria-controls', list.id);
    const optionId = (line: string) => `${ids}-${line.replace(/[^\w-]/g, '_')}`;
    // a model not in the atlas: planned from its own history, read in a dialog; the dropdown's foot, always in view
    const own = document.createElement('button');
    own.type = 'button';
    own.className = 'pmx-line-own';
    own.textContent = '+ Add your own model: plan it from its history…';
    own.title = 'A model the atlas has never seen, planned from its implantation pilot, the provider\'s growth data, ' +
      'or your earlier studies.';
    const panel = div([groups.root, search, list, own], 'pmx-line-panel');

    const visible = () => modelsFound(this._models.peek() ?? {}, cls.peek(), query.peek()).flatMap((g) => g.lines);
    const pick = (line: string) => {
      this.value = line;
      if (!line.startsWith(HISTORY_PREFIX))
        setUserSetting(LAST, line);
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

    this.build(() => {
      this.scope.listen(field, 'click', () => open.value = !open.peek());
      Tooltip.bind(note, () => FIT_NOTES[this.peek() ?? ''] ?? '');
      this.scope.listen(search, 'input', () => query.value = search.value);
      this.scope.listen(own, 'click', () => {
        open.value = false;
        historyDialog(pick);
      });
      this.scope.listen(panel, 'keydown', (ev) => {
        const e = ev as KeyboardEvent;
        if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
          e.preventDefault();
          move(e.key === 'ArrowDown' ? 1 : -1);
        } else if (e.key === 'Enter') {
          // what Enter picks is on screen: the marked line while the search still shows it, else the first shown
          const shown = visible();
          const line = shown.includes(active.peek() ?? '') ? active.peek()! : shown[0];
          if (!line)
            return;
          e.preventDefault();
          e.stopPropagation();
          pick(line);
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
        const own = historyModelName(v);
        // (none yet: a study read from a file names no model, and none is assumed)
        name.textContent = own ?? (v?.startsWith(HISTORY_PREFIX) ? 'Own history, not loaded' :
          v || 'Pick the tumor model');
        note.style.display = FIT_NOTES[v ?? ''] ? '' : 'none';
        note.setAttribute('aria-label', FIT_NOTES[v ?? ''] ?? '');
        kind.textContent = own ? 'Own history · PDX' :
          m ? `${m.tumour_type} · ${CLASSES.find((c) => c.cls === m.class)?.label ?? m.class}` : '';
        // the field's caption is its section's title on screen; a screen reader hears what the field is
        field.setAttribute('aria-label', `Tumor model: ${[name.textContent, kind.textContent].filter((x) => x)
          .join(', ')}`);
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
        if (at && models?.[at])
          search.setAttribute('aria-activedescendant', optionId(at));
        else
          search.removeAttribute('aria-activedescendant');
        if (!models) {
          list.replaceChildren(span('Loading the atlas’s models…', 'pmx-line-empty'));
          return;
        }
        const found = modelsFound(models, cls.value, query.value);
        // the models planned from their history in this session, to pick again after another (read each time the
        // list opens: one added meanwhile is there)
        void open.value;
        const words = query.value.toLowerCase().split(/\s+/).filter((w) => w);
        const mine = cls.value !== null && cls.value !== 'pdx' ? [] : historyKeys().filter((k) =>
          words.every((w) => (historyModelName(k) ?? '').toLowerCase().includes(w)));
        const ownItem = (key: string) => {
          const item = div([div([span(historyModelName(key) ?? key, 'pmx-line-item-name'),
            span('Own history · PDX', 'pmx-line-item-type')])], 'pmx-line-item');
          item.dataset.line = key;
          item.setAttribute('role', 'option');
          item.setAttribute('aria-selected', String(key === chosen));
          item.addEventListener('click', () => pick(key));
          return item;
        };
        const ownGroup = mine.length ? [Object.assign(span('Your models', 'pmx-line-group'), {role: 'presentation'}),
          ...mine.map(ownItem)] : [];
        // (a class's heading is no option: a screen reader reads it as text, not as a model to pick)
        const heading = (c: ModelClass) => {
          const e = span(CLASSES.find((x) => x.cls === c)!.label, 'pmx-line-group');
          e.setAttribute('role', 'presentation');
          return e;
        };
        list.replaceChildren(...ownGroup, ...found.length === 0 ? mine.length ? [] :
          [span('No model matches', 'pmx-line-empty')] :
          found.flatMap((g) => [heading(g.cls),
            ...g.lines.map((line) => {
              const m = models[line];
              const item = div([div([span(line, 'pmx-line-item-name'), span(m.tumour_type, 'pmx-line-item-type')]),
                span(`${m.papers} papers · ${m.animals} animals`, 'pmx-line-item-evidence')], 'pmx-line-item');
              item.dataset.line = line;
              item.id = optionId(line);
              item.setAttribute('role', 'option');
              item.setAttribute('aria-selected', String(line === chosen));
              item.classList.toggle('is-active', line === at);
              item.title = `${line}: ${m.tumour_type} (${m.organ}), ${m.host_class} host. Its growth prior comes ` +
                `from ${m.papers} papers, ${m.animals} animals.` + (FIT_NOTES[line] ? ` ${FIT_NOTES[line]}` : '');
              item.addEventListener('click', () => pick(line));
              return item;
            })]));
      });
    });
    return div([field, panel], 'pmx-line-input');
  }
}
