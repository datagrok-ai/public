import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {signal, computed, Scope, div, divV, divH, span, h3, link, iconButton, BoolInput, QNum}
  from '@datagrok-libraries/u2';
import {propertyForm, propertyEditor} from '@datagrok-libraries/u2/src/dg/index.js';
import {EDITOR_TYPES, enabledEditorTypes, registeredEditorTypes, setEditorTypes} from './editors';

const RUN = 'Local mode: `grok-core local stage U2Demo`, `grok-core local watch U2Demo`, ' +
  'open http://localhost:63343/login.html?mode=local. Dev stand: `npm run debug-u2demo`, ' +
  'then http://localhost:8082.';

const LOCAL = 'Two things local mode cannot do, both pre-existing and neither ours: the app card ' +
  'in Browse shows a red `error` chip (the platform asks the server for publication details that ' +
  'local mode has none of), and a reload cannot deep-link back — reopen `Apps / Dev / U2 Demo` ' +
  'and this tab comes back on its own.';

const HOVER = 'Number rows keep the platform\'s hover chrome, ported as it is: the slider and the ' +
  '− + clicker are invisible at rest and appear the moment the pointer enters the input. Each ' +
  'side reveals its OWN chrome and nothing else, so hover the two `Dose` cells one after the ' +
  'other — a bare cell means the pointer is on the other one, not that the chrome is missing. ' +
  'The revealed slider overlaps the row below by 9px — that is what the platform does, not a ' +
  'layout bug.';

const EDITOR = 'The leading `…` opens that row\'s own metadata in the context panel: an identity ' +
  'section, plus the fields its type actually has (bounds and slider flags for numbers, choices ' +
  'for strings). Every edit rebuilds BOTH editors of the whole grid from the new metadata — set ' +
  '`max` to 500 on `Dose` and watch both sides validate at 500, or switch `Name` from `string` ' +
  'to `int` and watch both editors change kind and drop the value.';

const NARROW = 'Two editor columns need about 700px between them. With Browse and the context ' +
  'panel both open the pane is narrower than that, and the grid — not the page — scrolls ' +
  'sideways: the prose stays put, and the row you are editing scrolls back with a drag of the ' +
  'bar under the matrix. Close the panel (or Browse) and both columns fit again.';

const CHEM_NOTE ='Structure (semType `Molecule`) needs the Chem package — both sides route to ' +
  'its sketcher, so the row is left out here. Open the app on a stand that has Chem published.';

const FILE_NOTE = 'Protocol (`file`) is degraded in local mode: dapi answers every existence ' +
  'check with nothing, so both editors settle on "does not exist". On a stand the row picks real ' +
  'files, with the platform\'s file browser behind both.';

const BIGINT_NOTE = 'Compound id (`bigint`) starts empty on purpose: a JS `bigint` in the target ' +
  'crashes the platform\'s own binding — `BigIntInput.value` calls `toString()` on it as if it ' +
  'were a Dart BigInt (`big_int_input.dart:29`). Type into either editor and both sides follow.';

const TAGS_NOTE = 'One row of the matrix is missing on purpose: `list` + inputType `Tags`. The ' +
  'platform reaches its `TagsInput` through the factory but binds only the value — nothing ' +
  'carries `choices` into `suggestionItems`, and `allowNew` stays false — so that editor can ' +
  'only delete the tag it was given, never add one (`tags_input.dart:23,374-386`, ' +
  '`input_base.dart:211-239`). u2 renders the row fine; the A/B would compare against a dead ' +
  'cell, so the row waits for a platform fix.';

/** The convergence matrix (plan.md §"Convergence page — the full input matrix"), in the order the
 * table numbers its rows: every row is one property-options record, rendered by the platform on
 * the left and by u2 on the right. */
function matrix(): DG.IProperty[] {
  return [
    {name: 'name', type: DG.TYPE.STRING, friendlyName: 'Name', nullable: false,
      description: 'Required on both sides'},
    {name: 'stage', type: DG.TYPE.STRING, friendlyName: 'Stage',
      choices: ['Discovery', 'Phase I', 'Phase II']},
    {name: 'route', type: DG.TYPE.STRING, friendlyName: 'Route', inputType: 'Radio',
      choices: ['Oral', 'Intravenous', 'Topical']},
    {name: 'notes', type: DG.TYPE.STRING, friendlyName: 'Notes', editor: 'textarea'},
    {name: 'apiKey', type: DG.TYPE.STRING, friendlyName: 'API key', editor: 'password'},
    {name: 'color', type: DG.TYPE.STRING, friendlyName: 'Color', inputType: 'Color'},
    {name: 'font', type: DG.TYPE.STRING, friendlyName: 'Font', inputType: 'Font'},
    {name: 'image', type: DG.TYPE.STRING, friendlyName: 'Image', inputType: 'Image'},
    {name: 'smiles', type: DG.TYPE.STRING, friendlyName: 'Structure', semType: 'Molecule'},
    {name: 'replicates', type: DG.TYPE.INT, friendlyName: 'Replicates', min: 1, max: 10},
    {name: 'cycles', type: DG.TYPE.INT, friendlyName: 'Cycles', min: 0, max: 50, showSlider: true},
    {name: 'dose', type: DG.TYPE.FLOAT, friendlyName: 'Dose', units: 'mg', min: 0, max: 1000,
      format: '#0.0', showSlider: true},
    {name: 'progress', type: DG.TYPE.NUM, friendlyName: 'Progress', editor: 'slider', min: 0,
      max: 100, step: 5},
    {name: 'compoundId', type: DG.TYPE.BIG_INT, friendlyName: 'Compound id'},
    {name: 'ic50', type: DG.TYPE.QNUM, friendlyName: 'IC50'},
    {name: 'active', type: DG.TYPE.BOOL, friendlyName: 'Active'},
    {name: 'archived', type: DG.TYPE.BOOL, friendlyName: 'Archived', editor: 'switch'},
    {name: 'started', type: DG.TYPE.DATE_TIME, friendlyName: 'Started'},
    {name: 'targets', type: DG.TYPE.LIST, friendlyName: 'Targets',
      choices: ['GPCR', 'Kinase', 'Ion channel', 'Protease']},
    {name: 'aliases', type: DG.TYPE.LIST, friendlyName: 'Aliases'},
    {name: 'attributes', type: DG.TYPE.MAP, friendlyName: 'Attributes'},
    {name: 'protocol', type: DG.TYPE.FILE, friendlyName: 'Protocol'},
  ];
}

/** One seed value per row, so both sides start on the same non-empty state. `dayjs`, not `Date`:
 * the platform reads a datetime property through its own converters, and a plain `Date` reaches
 * the Dart date editor as an unparseable object. */
function seed(): Record<string, any> {
  return {
    name: 'Aspirin',
    stage: 'Phase I',
    route: 'Oral',
    notes: 'Dissolved in saline, stored at 4 °C.',
    apiKey: 'sk-demo-0001',
    color: '#40a8e0',
    font: 'normal normal 12px "Roboto"',
    image: 'data:image/svg+xml;utf8,' + encodeURIComponent(
      '<svg xmlns="http://www.w3.org/2000/svg" width="64" height="40">' +
      '<rect width="64" height="40" rx="4" fill="#40a8e0"/></svg>'),
    smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
    replicates: 3,
    cycles: 12,
    dose: 250,
    progress: 40,
    // seeded empty on purpose: a JS `bigint` in the target crashes the platform's own binding —
    // `BigIntInput.value` calls `toString()` on it as if it were a Dart BigInt
    // (`big_int_input.dart:29`), and the direct property read skips js-api's marshalling
    compoundId: null,
    ic50: QNum.less(5.2),
    active: true,
    archived: false,
    started: dayjs('2026-01-15'),
    targets: ['Kinase'],
    aliases: ['ASA', 'acetylsalicylic acid'],
    attributes: {vendor: 'Sigma', lot: 'A-1234'},
    protocol: null,
  };
}

/** What a row's value resets to when its type is edited into another one. */
function defaultValue(type: string | undefined): any {
  switch (type) {
    case DG.TYPE.STRING:
      return '';
    case DG.TYPE.BOOL:
      return false;
    case DG.TYPE.LIST:
      return [];
    case DG.TYPE.MAP:
      return {};
    default:
      return null;
  }
}

/** Echo suppression for the cross-sync: a value that already matches is never written back.
 * Datetimes cross the bridge as `dayjs` one way and `Date` the other, so they compare as ms;
 * lists and maps are rebuilt on every read, so they compare by content. */
function same(a: any, b: any): boolean {
  if (a == null || b == null)
    return (a == null) === (b == null);
  const ta = time(a);
  const tb = time(b);
  if (ta !== undefined && tb !== undefined)
    return ta === tb;
  if (typeof a === 'number' && typeof b === 'number')
    return Object.is(a, b);
  if (Array.isArray(a) && Array.isArray(b))
    return JSON.stringify(a) === JSON.stringify(b);
  if (isRecord(a) && isRecord(b))
    return JSON.stringify(a) === JSON.stringify(b);
  return a === b;
}

/** A plain object — a `FileInfo` or any other class instance is compared by identity instead. */
function isRecord(x: any): boolean {
  return typeof x === 'object' && Object.getPrototypeOf(x) === Object.prototype;
}

function time(x: any): number | undefined {
  return x instanceof Date || typeof x?.toDate === 'function' ? +x.valueOf() : undefined;
}

/** Page prose with `backticked` spans rendered as code rather than shown as backticks. */
function prose(text: string): HTMLElement {
  const line = span('', 'u2demo-hint');
  const parts = text.split('`');
  for (let i = 0; i < parts.length; i++)
    line.append(i % 2 === 0 ? parts[i] : span(parts[i], 'u2demo-code'));
  return line;
}

/** The page a `PropRow` belongs to — it owns the one editor the context panel shows. */
interface PropPage {
  editorFor(name: string): HTMLElement;
}

/** What the '…' column makes current. */
export class PropRow {
  constructor(readonly name: string, readonly page: PropPage) {}
}

/** Pushes a custom editor into the context panel the way DevTools' signature editor does:
 * an ObjectHandler over a holder class, `grok.shell.o` as the door
 * (`DevTools/src/signature-editor/function-signature-editor.ts:35-62,302`). */
class PropRowHandler extends DG.ObjectHandler<PropRow> {
  get type(): string {
    return 'u2demo-prop-editor';
  }

  isApplicable(x: any): boolean {
    return x instanceof PropRow;
  }

  renderProperties(x: PropRow): HTMLElement {
    return x.page.editorFor(x.name);
  }
}

let handlerRegistered = false;

/** Registered once per session: the platform keeps every handler it is given, and reopening the
 * app would otherwise stack duplicates. */
export function registerPropRowHandler(): void {
  if (handlerRegistered)
    return;
  handlerRegistered = true;
  DG.ObjectHandler.register(new PropRowHandler());
}

export function convergencePage(): HTMLElement {
  const scope = Scope.ambient!;
  const hasChem = DG.Func.find({package: 'Chem', name: 'drawMolecule'}).length > 0;
  const target = seed();
  let records: Record<string, DG.IProperty> = {};
  for (const record of matrix())
    records[record.name!] = record;

  const toDart = signal(0);
  const toU2 = signal(0);
  // the platform delivers its own onChanged from a microtask, so the guard against counting our
  // own write-through as a Dart edit has to outlive the synchronous application below
  let applying = 0;
  let edited: string | null = null;
  let gridScope: Scope | undefined;
  let scheduled = false;

  const grid = div([], 'u2demo-ab');
  scope.own(() => gridScope?.dispose());

  const header = h3('');
  const editor = propertyEditor(records[Object.keys(records)[0]], {
    onChanged: (record) => {
      if (edited == null)
        return;
      edited = apply(edited, record);
      header.textContent = `Property: ${edited}`;
      markCurrent();
      schedule();
    },
    // what `apply` refuses, said on the field instead of dropped in silence
    validators: {name: (x) => x !== edited && records[x] != null ? 'Name is already taken' : null},
  });
  const panel = divV([header, editor.root]);

  const owner: PropPage = {
    editorFor: (name) => {
      edited = name;
      header.textContent = `Property: ${name}`;
      editor.setTarget(records[name]);
      markCurrent();
      return panel;
    },
  };

  /** The row the context panel is editing, so the '…' that opened it is not the only clue. */
  function markCurrent(): void {
    const cells = grid.querySelectorAll<HTMLElement>('[data-row]');
    for (let i = 0; i < cells.length; i++)
      cells[i].classList.toggle('u2demo-ab-current', cells[i].dataset.row === edited);
  }

  /** The editor reports the whole record, `name` included: a rename re-keys the record in place
   * and carries the value over, so the row keeps its position. A name that is empty or already
   * taken is not applied; the editor keeps what was typed and says why on the field. */
  function apply(key: string, record: DG.IProperty): string {
    const name = record.name ?? '';
    if (name !== key && (name === '' || records[name] != null))
      return key;
    const before = records[key];
    const next: Record<string, DG.IProperty> = {};
    for (const [k, value] of Object.entries(records))
      next[k === key ? name : k] = k === key ? record : value;
    records = next;
    if (name !== key) {
      target[name] = target[key];
      delete target[key];
    }
    if (record.type !== before.type)
      target[name] = defaultValue(record.type);
    return name;
  }

  function schedule(): void {
    if (scheduled)
      return;
    scheduled = true;
    queueMicrotask(() => {
      scheduled = false;
      build();
    });
  }

  /** What actually scrolls around the grid — the tab strip's content box, not the page div. */
  function scroller(): HTMLElement | null {
    for (let el = grid.parentElement; el != null; el = el.parentElement)
      if (el.scrollHeight > el.clientHeight + 1)
        return el;
    return null;
  }

  function build(): void {
    const scrolled = scroller();
    const top = scrolled?.scrollTop ?? 0;
    gridScope?.dispose();
    const generation = new Scope();
    gridScope = generation;
    grid.textContent = '';
    grid.append(span(''), h3('Dart — ui.input.forProperty'), h3('u2 — propertyForm'));

    Scope.runWith(generation, () => {
      const shown = Object.values(records).filter((r) => hasChem || r.semType !== 'Molecule');
      const props = shown.map((r) => DG.Property.fromOptions(r));
      const dart = props.map((p) => ui.input.forProperty(p, target));
      const form = propertyForm(props, target, {
        onChanged: () => {
          toDart.value++;
          applying++;
          try {
            for (let i = 0; i < props.length; i++) {
              const value = props[i].get(target);
              if (same(dart[i].value, value))
                continue;
              // toDart() converts dayjs but not a plain Date, which is what the u2 date editor holds
              dart[i].value = props[i].propertyType === DG.TYPE.DATE_TIME && value != null ?
                dayjs(value) : value;
            }
          } finally {
            queueMicrotask(() => applying--);
          }
        },
      });
      const refresh = () => {
        if (applying > 0)
          return;
        toU2.value++;
        form.refresh();
      };

      let i = 0;
      for (const record of Object.values(records)) {
        if (!hasChem && record.semType === 'Molecule') {
          grid.append(note(CHEM_NOTE));
          continue;
        }
        const name = record.name!;
        // the panel is closed by default in local mode, and shell.o alone would not open it
        const cells = [iconButton('ellipsis-h', () => {
          grok.shell.windows.showContextPanel = true;
          grok.shell.o = new PropRow(name, owner);
        }, {tooltip: 'Edit this property\'s metadata'}), dart[i].root, form.input(name)!.root];
        for (const cell of cells)
          cell.dataset.row = name;
        grid.append(...cells);
        if (record.type === DG.TYPE.FILE)
          grid.append(note(FILE_NOTE));
        if (record.type === DG.TYPE.BIG_INT)
          grid.append(note(BIGINT_NOTE));
        // where the missing `list` + `Tags` row would sit, right after the other list rows
        if (name === 'aliases')
          grid.append(note(TAGS_NOTE));
        i++;
      }
      markCurrent();

      for (const input of dart) {
        const sub = input.onChanged.subscribe(refresh);
        generation.own(() => sub.unsubscribe());
      }
    });
    if (scrolled)
      scrolled.scrollTop = top;
  }

  const enabled = new Set(enabledEditorTypes());
  const active = new Set(registeredEditorTypes());
  const toggles = EDITOR_TYPES.map((type) => new BoolInput({
    label: active.has(type) ? `${type} (active)` : type,
    value: enabled.has(type),
    onChanged: (on) => {
      if (on)
        enabled.add(type);
      else
        enabled.delete(type);
      setEditorTypes([...enabled]);
      grok.shell.info(on ? `u2 ${type} editor will register on next reload` :
        `${type} goes back to the Dart editor on next reload`);
    },
  }));

  const page = divV([
    prose('One property set, one target object, two generators: the platform\'s ' +
      '`InputBase.forProperty` on the left, u2\'s `propertyForm` on the right — every input type ' +
      'the platform routes to, one per row. Editing either side writes the same object and the ' +
      'other side follows. The pickers that need an open table — columns, columns map, ' +
      'aggregations, tables — are the same convergence work and live on the `Files & columns` tab.'),
    prose(HOVER),
    prose(EDITOR),
    prose('Invalid states diverge deliberately: the platform reddens the caption and the ' +
      'underline and says no more, u2 spells the reason out beside the field. Clear `Name` or ' +
      'type `99` into `Replicates` to see both.'),
    prose('Display formats diverge too. Both sides format on a programmatic write and leave what ' +
      'you are typing alone, so the side you edit shows `500` while the other shows `500.0`. ' +
      'Dates differ by design: the platform prints them in the locale format, u2 in ISO. So do ' +
      'qualified numbers: `IC50` reads `<5.20` on the platform, which pads to its default ' +
      'precision, and `<5.2` in u2, which prints the value as it stands.'),
    prose(NARROW),
    prose(RUN),
    prose(LOCAL),
    grid,
    divH([span('sync = '), span(computed(() =>
      `u2 → Dart ${toDart.value} · Dart → u2 ${toU2.value}`))], 'u2demo-status'),
    h3('Dual run: u2 inputs as platform value editors'),
    prose('Each type registers a `role: valueEditor` func, which the platform consults before any ' +
      'type default — so the LEFT column, still built by the platform, starts rendering u2 ' +
      'editors for that type. Registration happens once per session: types marked `(active)` ' +
      'registered when this session started, the rest apply on the next reload — including an ' +
      'unchecked one going back to the Dart editor. Enable `string` last — it also captures ' +
      'pattern and password properties, which u2 does not reimplement.'),
    divH([...toggles, link('Reload now', () => window.location.reload())], 'u2demo-toggles'),
  ], 'u2demo-page u2demo-wide');

  build();
  return page;
}

/** A full-width row in the A/B grid carrying a degradation note for the row above it. */
function note(text: string): HTMLElement {
  const el = prose(text);
  el.classList.add('u2demo-ab-note');
  return el;
}
