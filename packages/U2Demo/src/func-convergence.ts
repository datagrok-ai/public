import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {signal, computed, Scope, div, divV, divH, span, h3} from '@datagrok-libraries/u2';
import {funcForm, FuncCallForm, ObjectForm, IProperty} from '@datagrok-libraries/u2/src/dg/index.js';
import {kindOf} from '@datagrok-libraries/u2/src/dg/forms/object-form.js';
import type {Kind} from '@datagrok-libraries/u2/src/dg/forms/object-form.js';

const RUN = 'Local mode: `grok-core local stage U2Demo`, `grok-core local watch U2Demo`, ' +
  'open http://localhost:63343/login.html?mode=local. Dev stand: `npm run debug-u2demo`, ' +
  'then http://localhost:8082.';

const INTRO = 'One FuncCall, two editors: the platform\'s `DG.InputForm.forFuncCall` on the left, ' +
  'u2\'s `funcForm` on the right — the W1 scalar matrix (static choices with a nullable empty ' +
  'option, radio, textarea, color, slider, bounded numbers with units, plain bools, datetime), ' +
  'generated from one client-registered function decorated after registration. Editing either ' +
  'side writes the same FuncCall and the other side follows; the `inputs` line below is what the ' +
  'call itself holds. Each cross-edit bumps exactly one sync counter by exactly one — an echo ' +
  'never counts, so any extra increment is ping-pong.';

const DEFAULTS = 'Literal numeric defaults ride the js-api `defaultValue` setter display-only ' +
  'and, today, u2-only: `Replicates` and `Dose Level` open showing 3 and 250 on the u2 side ' +
  'while the platform form leaves them empty — the platform func form reads defaults only from ' +
  '`options[\'default\']` (`func_param_editor.dart:558`), a GrokScript command it evaluates and ' +
  'writes into the call. Either way `inputs` reads null until you edit.';

const STYLE = '`Color` sits under a `Style` header neither side asked for: with no category set, ' +
  'the platform derives one from the NAME (`property.dart:340` — anything ending in `color`), ' +
  'and both sides read the same property, so they group identically.';

const VALID = 'Invalid states diverge deliberately: the platform reddens the caption and the ' +
  'underline and says no more, u2 spells the reason out beside the field. Type `99` into ' +
  '`Replicates` to see both — out of range is a validation error on both sides, never a clamp. ' +
  'The invalid value still lands in the call — parity: validation gates the run, not the write.';

const RADIO = 'Dart-side radio edits never reach the call (`RadioInput` fires `onInput` only; ' +
  'the func form writes on `onChanged`) — edit `route` on the u2 side; Dart follows.';

const SLIDERS = 'At rest neither slider was written: Dart shows the browser range default (50), ' +
  'u2 the range floor (0) — the call holds `progress: null` until a drag.';

const DATETIME = 'Datetime crosses the bridge as `dayjs`; u2 holds a `Date` and converts at the ' +
  'write boundary — cross-edits settle without ping-pong (watch the counters stay put). Display ' +
  'formats differ by design: the platform prints the date in the locale format, u2 in ISO.';

const REQUIRED = 'The Dart pattern-parameter exemption from required validation is out of W1 ' +
  'scope; the required rule here is a plain not-empty check.';

const FRIENDLY = 'The `friendlyName` decoration is u2-only for now: the platform func form ' +
  'captions the input from it, then looks the input back up by the caption derived from the NAME ' +
  '(`func_param_editor.dart:834`) — any friendlyName-decorated scalar param nulls that lookup ' +
  'and the whole form throws (`addValidators` on null). This bench names the param `doseLevel` ' +
  'instead, which camelCases to `Dose Level` on both sides; the u2 form keys its lookups by ' +
  'name to avoid exactly this wart.';

const STRING_DEFAULT = 'A string `defaultValue` set through the js-api lands in ' +
  '`options[\'default\']`, which the platform evaluates as a GrokScript command — `"Dissolved ' +
  'in saline."` balloons `Unable to calculate default value`. Numeric defaults stay literal and ' +
  'display-only on both sides.';

const GATE = 'A second function with one required parameter: `isValid` is read off both forms ' +
  'after every edit. Type into `batchId` and both turn true; clear it and both turn false. ' +
  'Before the first edit they disagree on purpose to observe: the platform answers valid until ' +
  'something makes it validate, u2 validates eagerly.';

const DYN_INTRO = 'One shared FuncCall again, but every source is a function: `country` lists ' +
  '`fceW2Countries()`, `city` re-asks `fceW2Cities(country)` after a country edit on EITHER ' +
  'side, `model` comes off a DataFrame provider whose picked row fills `mpg`/`cyl` ' +
  '(`propagateChoice: \'all\'`), `site` suggests from the typed text, `cohort` opens with its ' +
  'computed default already written into the call, `seed`\'s default command is broken on ' +
  'purpose, and `regions` runs the country provider through a multi-choice. `country` also ' +
  'carries `descriptions` — hover the input for the selected item\'s.';

const DEP_REFRESH = 'Change `country` on either side and both `city` dropdowns re-ask the ' +
  'provider. u2 debounces every dependent re-eval at a unified 200 ms — the Dart editor ' +
  'debounces list params only and re-runs string choices on every dependency edit — so a ' +
  'propagate burst coalesces into one re-eval per source.';

const PRUNE = 'Pick a city, then change `country`: the refreshed items no longer offer it, and ' +
  'u2 prunes the stale value INTO the call — `city` reads null in `inputs` below. The Dart ' +
  'select goes blank but the call keeps holding the invisible value (divergence #9).';

const PROPAGATE = 'Pick a model on either side and that row\'s `mpg`/`cyl` land in the call: ' +
  'u2 writes them through the sibling inputs, so their validation runs and a dependent source ' +
  'downstream would re-fire once. The provider columns are ints on purpose — float32 noise is ' +
  'what a double column would add.';

const SUGGEST = 'Both sides now ask `fceW2Suggest` with the TYPED text (the W1 Dart-side fix; ' +
  'the pre-fix editor evaluated the stored value). u2 renders the popup in its SuggestionList ' +
  'vocabulary — loading, empty and error-with-Retry rows the platform popup does not have.';

const DEFAULT_OK = '`options[\'default\'] = \'2 + 2\'` is a GrokScript command evaluated and ' +
  'WRITTEN into the call at open (R6) — `cohort: 4` sits in `inputs` below before any edit. ' +
  'Defaults run once per bind and never re-run on dependency edits: a default command cannot ' +
  'reference sibling params.';

const DEFAULT_BROKEN = '`seed`\'s default command names a variable that does not exist: the ' +
  'platform balloons `Unable to calculate default value` at open and moves on; u2 keeps the ' +
  'failure on the field — the message inline, with a Retry that re-runs just that one eval.';

const LIST_DYN = 'The same provider through a multi-choice on both sides. Under ' +
  '`skipDefaultInit` the Dart guard (`func_param_editor.dart:622`) skips the ITEMS write too ' +
  'and leaves the dropdown permanently empty; u2 applies freshly loaded items either way — ' +
  'items are not defaults (divergence #8).';

const W3_INTRO = 'One shared call over a `DG.Script.create` fixture — annotation-true metadata ' +
  '(`column num {type: numerical; caption: Numeric column}` — the caption reaches both labels — ' +
  '`column cat {type: categorical}`, `column_list metrics ' +
  '{type: numerical}`, all implicitly linked to `df`), no server round-trip. Two fixture tables ' +
  'are opened for it: `fceW3Demog` (age, height, weight, sex) and `fceW3Alt` (weight, label). ' +
  'Both forms auto-fill the SAME table into the call at open, and `num`/`cat` are auto-picked ' +
  'by name similarity — real writes, visible in `inputs` below before any edit. Column lists ' +
  'are never auto-picked. The call holds OBJECTS (a raw name write comes back as a pending ' +
  'resolver call), the inputs hold names. One deliberate exception to the one-tick counter rule: ' +
  'a table pick also re-defaults every dependent through its input, so the sync counter counts ' +
  'the pick plus each re-defaulted column — a 3 here is expected, not ping-pong.';

const W3_ORDER = 'Declaration order (divergence #10): u2 renders `df, num, cat, metrics` as ' +
  'declared; the raw Dart form clusters dependents at the table\'s slot (`Df, Metrics, Num, ' +
  'Cat`). This grid re-parents both sides into declaration-order rows, so the difference shows ' +
  'only outside the bench.';

const W3_AUTOPICK = 'The auto-pick runs regardless of `skipDefaultInit` — the Dart guard covers ' +
  'only `options[\'default\']` (`func_param_editor.dart:558`) and u2 reproduces that. Picks ' +
  'follow the Dart distance algorithms (Levenshtein for 1-character param names, Jaro-Winkler ' +
  'otherwise); a semType-carrying param takes the molecule-named or first passing column instead.';

const W3_NOMATCH = 'Switch to a table with no categorical column and `cat` re-defaults to NULL ' +
  'written into the call, the u2 select left empty — Dart itself prunes to the call here, so ' +
  'this is exact parity, not a divergence.';

const W3_EXTERNAL = 'u2 subscribes the table PARAM, one subscription for all dependents ' +
  '(divergence #12): an external `setParamValue` of another table retargets the column pickers ' +
  'too. Dart wires the table INPUT and stacks one listener per dependent, so external table ' +
  'writes never rewire its combos.';

const W3_CLOSED = 'Close the table the call holds and u2 prunes the vanished value INTO the ' +
  'call — `df` and the dependent columns read null below (divergence #11). The Dart form keeps ' +
  'the closed frame in the call while its select silently displays a remaining table it does ' +
  'not hold.';

const W3_LIST = 'Pick columns through the u2 popup and a `DG.Column[]` goes into the call, ' +
  'reading back as a `ColumnList` — `(2) …` on the Dart side, the names in `inputs` below. A ' +
  'table switch clears the selection to `[]`; a pre-set value survives form construction. The ' +
  'popup itself is now the platform ColumnsGrid in checkbox mode with dialog semantics: checks ' +
  'buffer in the grid until OK commits ONE write; Cancel, Esc and an outside click discard them.';

const W3_UNASSOC = 'Not on this grid by construction: a column param whose table param cannot ' +
  'be resolved is LISTED under `unsupported` by u2 (divergence #13) — the Dart form silently ' +
  'renders nothing for it. u2 resolves the annotation link first (`parentTableParamName`), then ' +
  'the explicit `{table: …}` option, which also serves signature-registered funcs decorated at ' +
  'runtime; `editor: columnsMap` params stay unsupported — the shipped Dart form crashes on ' +
  'them (platform defect #9).';

const P_ICONS = 'The u2 table field now carries the platform rail: \'Open file\' plus the ' +
  'feature-detected \'Add file from Files\' and \'Query database\' icons, opening the REAL ' +
  'platform dialogs through `ui.pickTableFromFiles`/`ui.pickTableFromQuery` — the promise ' +
  'settles null on Cancel or Esc and the form stays untouched. One divergence (#16): a picked ' +
  'frame joins the WORKSPACE and becomes the value as a user edit — the Dart icons feed the ' +
  'input\'s private item list only, invisible to every other table input. The rail is ' +
  'per-instance: read-modify-write the `TableInput.actions` list.';

const P_COMBO = 'The column fields are now the grid combo: click (or type on) either one and ' +
  'the platform\'s own ColumnsGrid drops under the field — the same `d4-column-grid` the Dart ' +
  '`d4-column-selector` on the left opens, anchored in a u2 overlay, so both sides share one ' +
  'dropdown experience. Three deliberate divergences: the search box is visible from open ' +
  'instead of type-to-reveal (#15), hovering a row never previews the value into the call ' +
  '(#14 — in a form each preview would fire the dependent cascade), and the mouse wheel does ' +
  'not cycle the value over the closed field (#17) — the arrow keys still do. A popup pick is ' +
  'a user edit even when it confirms the current value, so it clears the `auto` badge.';

const DIVERGENCES = [
  '`editor:` hints (`textarea`, `password`, `switch`, `slider`) — u2 honors them, type-guarded ' +
  'the way the property form does; the platform func form evaluates any such hint as a nested ' +
  'function editor and throws (`Variable "textarea" not found`), so nothing can depend on them ' +
  'there. This bench therefore decorates via `inputType` only — and since no `Password` ' +
  'inputType exists, `secret` stays a plain text row on both sides for now.',
  'Suggestions: the platform editor used to evaluate `suggestions` against the parameter\'s ' +
  'stored value instead of the typed text. Fixed Dart-side (`evalParamSuggestions` carries the ' +
  'text) and consumed here: the u2 typeahead on `site` queries with what was typed and shows ' +
  'its own loading, empty and error-with-Retry rows.',
  'Lookups are name-keyed in u2: `getInput(\'doseLevel\')`, never `getInput(\'Dose Level\')` — ' +
  'the platform\'s caption-keyed lookups are not reproduced.',
  'One subscription per parameter, dropped on a source rebind — no listener accumulation, and ' +
  'no `needsValidation` memo: u2\'s validity is a live signal.',
];

const SIGNATURE = 'string fceDemo(string stage, string route, string notes, string secret, ' +
  'string color, int replicates, double doseLevel, int progress, bool active, bool archived, ' +
  'datetime started)';

let demoFunc: DG.Func | null = null;
let gateFunc: DG.Func | null = null;
let dynFunc: DG.Func | null = null;
let w3Func: DG.Func | null = null;

/** Opens a named fixture table when the workspace does not already hold it — idempotent, so a
 * reopened tab (or a user who closed the table) gets it back; `addTable` would dedupe the name
 * into a copy otherwise. */
export function ensureTable(name: string, csv: string): DG.DataFrame {
  const open = grok.shell.tables.find((t) => t.name === name);
  if (open != null)
    return open;
  const df = DG.DataFrame.fromCsv(csv);
  df.name = name;
  grok.shell.addTable(df);
  return df;
}

/** Registered once per session (the `registerPropRowHandler` pattern): the platform keeps every
 * function it is given, and reopening the app would otherwise stack duplicates. The signature
 * parser takes only `type name` pairs, so every decoration lands after registration, through the
 * js-api Property setters and the live options map. */
export function ensureFuncs(): void {
  if (demoFunc != null)
    return;
  demoFunc = grok.functions.register({
    signature: SIGNATURE,
    run: (...args: any[]) => args.map((x) => String(x ?? '')).join(' '),
  });
  const inputs = new Map(demoFunc.inputs.map((p) => [p.name, p] as [string, DG.Property]));
  const stage = inputs.get('stage')!;
  stage.choices = ['Discovery', 'Phase I', 'Phase II'];
  stage.nullable = true;
  const route = inputs.get('route')!;
  route.inputType = 'Radio';
  route.choices = ['Oral', 'Intravenous', 'Topical'];
  const notes = inputs.get('notes')!;
  notes.inputType = 'TextArea';
  notes.description = 'Free-form dosing notes';
  const color = inputs.get('color')!;
  color.inputType = 'Color';
  // signature-registered params come out nullable: false; without this the matrix rests under
  // a wall of "Value can't be empty" — required-ness is the gate function's row to demonstrate
  for (const name of ['route', 'notes', 'secret', 'color', 'started'])
    inputs.get(name)!.nullable = true;
  const replicates = inputs.get('replicates')!;
  replicates.min = 1;
  replicates.max = 10;
  replicates.defaultValue = 3;
  const dose = inputs.get('doseLevel')!;
  dose.min = 0;
  dose.max = 1000;
  dose.showSlider = true;
  dose.defaultValue = 250;
  dose.options['units'] = 'mg';
  const progress = inputs.get('progress')!;
  progress.inputType = 'Slider';
  progress.min = 0;
  progress.max = 100;
  progress.category = 'Advanced';
  const started = inputs.get('started')!;
  started.category = 'Advanced';
  started.description = 'When the study started';
  gateFunc = grok.functions.register({
    signature: 'bool fceGate(string batchId)',
    run: (batchId: string) => batchId != null && batchId !== '',
  });
  gateFunc.inputs[0].nullable = false;
  grok.functions.register({
    signature: 'List<String> fceW2Countries()',
    run: () => ['FR', 'DE', 'US'],
  });
  grok.functions.register({
    signature: 'List<String> fceW2Cities(String country)',
    run: (c: string) => c == null || c === '' ? [] : [`${c}-1`, `${c}-2`],
  });
  grok.functions.register({
    // int columns on purpose: a double column round-trips with float32 noise (22.8 → 22.7999…)
    signature: 'dataframe fceW2CarsDf()',
    run: () => DG.DataFrame.fromCsv('model,mpg,cyl\nMazda RX4,21,6\nDatsun 710,22,4\nHornet 4 Drive,18,8'),
  });
  grok.functions.register({
    signature: 'List<String> fceW2Suggest(String s)',
    run: (s: string) => [`${s}-1`, `${s}-2`],
  });
  dynFunc = grok.functions.register({
    signature: 'string fceW2Dyn(string country, string city, string model, int mpg, int cyl, ' +
      'string site, int cohort, int seed, list regions)',
    run: (...args: any[]) => args.map((x) => String(x ?? '')).join(' '),
  });
  const dyn = new Map(dynFunc.inputs.map((p) => [p.name, p] as [string, DG.Property]));
  const country = dyn.get('country')!;
  country.options['choices'] = 'fceW2Countries()';
  country.options['descriptions'] = {FR: 'France', DE: 'Germany', US: 'United States'};
  dyn.get('city')!.options['choices'] = 'fceW2Cities';
  const model = dyn.get('model')!;
  model.options['choices'] = 'fceW2CarsDf()';
  model.options['propagateChoice'] = 'all';
  dyn.get('site')!.options['suggestions'] = 'fceW2Suggest';
  dyn.get('cohort')!.options['default'] = '2 + 2';
  dyn.get('seed')!.options['default'] = 'nosuchvar + 1';
  dyn.get('regions')!.options['choices'] = 'fceW2Countries()';
  for (const p of dynFunc.inputs)
    p.nullable = true;
  // Script.create, not a signature registration: the annotation parser is what sets the implicit
  // table link and the column filters — the Dart-parity path both forms must share
  w3Func = DG.Script.create([
    '//name: fceW3Bench',
    '//language: javascript',
    '//input: dataframe df',
    '//input: column num {type: numerical; caption: Numeric column}',
    '//input: column cat {type: categorical}',
    '//input: column_list metrics {type: numerical}',
    'return;',
  ].join('\n'));
}

export function fmt(v: any): string {
  if (v == null)
    return 'null';
  if (Array.isArray(v))
    return `[${v.join(', ')}]`;
  if (v instanceof DG.DataFrame)
    return `DataFrame '${v.name}', ${v.rowCount} rows × ${v.columns.length} columns`;
  if (v instanceof DG.Column)
    return `Column '${v.name}'`;
  if (typeof v?.names === 'function')
    return `[${v.names().join(', ')}]`;
  if (typeof v?.toDate === 'function')
    return v.toDate().toISOString();
  return typeof v === 'string' ? `'${v}'` : String(v);
}

function inputsText(call: DG.FuncCall): string {
  const parts: string[] = [];
  for (const p of call.inputParams.values())
    parts.push(`${p.name}: ${fmt(p.value)}`);
  return `{${parts.join(', ')}}`;
}

/** Page prose with `backticked` spans rendered as code rather than shown as backticks. */
export function prose(text: string): HTMLElement {
  const line = span('', 'u2demo-hint');
  const parts = text.split('`');
  for (let i = 0; i < parts.length; i++)
    line.append(i % 2 === 0 ? parts[i] : span(parts[i], 'u2demo-code'));
  return line;
}

export function funcConvergencePage(): HTMLElement {
  const scope = Scope.ambient!;
  ensureFuncs();
  // before any form builds: both sides' table auto-fill reads the open-tables list
  ensureTable('fceW3Demog', 'age,height,weight,sex\n25,170,70,M\n31,166,64,F\n44,182,80,M');
  ensureTable('fceW3Alt', 'weight,label\n70,x\n80,y');
  const call = demoFunc!.prepare();
  const gateCall = gateFunc!.prepare();
  // country preset so the dependent city source has something to ask with at open
  const dynCall = dynFunc!.prepare({country: 'FR'});
  const w3Call = w3Func!.prepare();

  const toDart = signal(0);
  const toU2 = signal(0);
  const tick = signal(0);
  const toDartDyn = signal(0);
  const toU2Dyn = signal(0);
  const tickDyn = signal(0);
  const toDartW3 = signal(0);
  const toU2W3 = signal(0);
  const tickW3 = signal(0);
  const fields = signal('…');
  const gateStatus = signal('…');

  const grid = div([], 'u2demo-ab u2demo-ab-func');
  const gateGrid = div([], 'u2demo-ab u2demo-ab-func');
  const dynGrid = div([], 'u2demo-ab u2demo-ab-func');
  const w3Grid = div([], 'u2demo-ab u2demo-ab-func');
  const plain = div([], 'u2demo-standalone');

  /** One A/B section: both editors over one shared call, re-parented into a grid row per param
   * with its notes, and a per-section honest counter pair. */
  async function buildAb(abCall: DG.FuncCall, abGrid: HTMLElement,
    toDartCount: {value: number}, toU2Count: {value: number}, tickCount: {value: number},
    rowNotes: Map<string, string[]>, catNotes?: Map<string, string>):
    Promise<{dart: DG.InputForm, form: FuncCallForm}> {
    const dartForm = await DG.InputForm.forFuncCall(abCall, {twoWayBinding: true});

    // the honest Dart → u2 counter: subscribed BEFORE funcForm, so each handler still sees the
    // field value the form's own (later-subscribed) refresh is about to replace; the comparison
    // is the refresh's own coerce + same, so it counts exactly the refreshes that change a field
    // — and a u2-originated echo, whose field already holds the value, never counts
    let form: FuncCallForm | undefined;
    const reads = new Map<string, (v: any) => any>();
    for (const p of abCall.inputParams.values()) {
      const prop = p.property as unknown as IProperty & {options?: Record<string, any>};
      const type = prop.propertyType ?? prop.type;
      // the W3 routes hold OBJECTS in the call and names in the inputs — compare names, the way
      // the form's own fromParam converters do (a pending resolver readback reads as null)
      if (type === 'dataframe' || type === 'column')
        reads.set(p.name, (v) =>
          v != null && typeof v.name === 'string' && !('func' in v) ? v.name : null);
      else if (type === 'column_list')
        reads.set(p.name, (v) => v == null ? [] : typeof v.names === 'function' ? v.names() :
          Array.isArray(v) ?
            v.map((c: any) => c?.name).filter((n: any) => typeof n === 'string') : []);
      else {
        // a dynamic-choices param has no typed prop.choices, so kindOf would answer 'string',
        // whose coerce turns null into '' — mirror the form's own routing instead
        const kind: Kind = type === 'string' && prop.options?.['choices'] != null ?
          'choice' : kindOf(prop, true);
        reads.set(p.name, (v) => ObjectForm.coerce(kind, v));
      }
      const sub = p.onChanged.subscribe(() => {
        tickCount.value++;
        const input = form?.getInput(p.name);
        if (input == null)
          return;
        const value = reads.get(p.name)!(p.value);
        if (!ObjectForm.same(value, input.value.peek()))
          toU2Count.value++;
      });
      scope.own(() => sub.unsubscribe());
    }

    form = Scope.runWith(scope, () => funcForm(abCall, {
      twoWayBinding: true,
      onInputChanged: () => toDartCount.value++,
    }));

    // one grid row per param, both sides' inputs re-parented into it (the convergence-page
    // pattern), grouped the way both forms group: categories in first-appearance order
    const categories: {name: string, params: DG.FuncCallParam[]}[] = [];
    for (const p of abCall.inputParams.values()) {
      const name = (p.property.category as string | null | undefined) ?? 'Misc';
      let category = categories.find((c) => c.name === name);
      if (category === undefined) {
        category = {name, params: []};
        categories.push(category);
      }
      category.params.push(p);
    }
    const dartByName = new Map<string, DG.InputBase>();
    for (const input of dartForm.inputs) {
      const prop = (input as any).property;
      if (prop?.name != null)
        dartByName.set(prop.name, input);
    }
    abGrid.append(span(''), h3('Dart — DG.InputForm.forFuncCall'), h3('u2 — funcForm'));
    for (const category of categories) {
      // the default group gets no invented header — only real categories announce themselves
      if (category.name !== 'Misc') {
        abGrid.append(div([category.name], 'u2demo-ab-cat'));
        const catNote = catNotes?.get(category.name);
        if (catNote != null)
          abGrid.append(note(catNote));
      }
      for (const p of category.params) {
        const u2Input = form.getInput(p.name);
        const dartInput = dartByName.get(p.name);
        if (u2Input == null || dartInput == null)
          continue;
        const cells = [span(p.name, 'u2demo-code'), dartInput.root, u2Input.root];
        for (const cell of cells)
          cell.dataset.row = p.name;
        abGrid.append(...cells);
        for (const text of rowNotes.get(p.name) ?? [])
          abGrid.append(note(text));
      }
    }
    return {dart: dartForm, form};
  }

  async function build(): Promise<void> {
    const {dart: dartForm, form} = await buildAb(call, grid, toDart, toU2, tick,
      new Map<string, string[]>([
        ['route', [RADIO]],
        ['notes', [STRING_DEFAULT]],
        ['replicates', [VALID]],
        ['doseLevel', [DEFAULTS, FRIENDLY]],
        ['progress', [SLIDERS]],
        ['started', [DATETIME]],
      ]), new Map([['Style', STYLE]]));
    fields.value = `Dart ${dartForm.inputs.length} · u2 ${form.inputs.length} · unsupported: ` +
      (form.unsupported.length === 0 ? '(none)' : form.unsupported.join(', '));

    const dartGate = await DG.InputForm.forFuncCall(gateCall, {twoWayBinding: true});
    const gateForm = Scope.runWith(scope, () => funcForm(gateCall, {
      twoWayBinding: true,
      onInputChanged: () => refreshGate(),
    }));

    function refreshGate(): void {
      // the platform validates from its own change handler; read both a microtask later
      queueMicrotask(() => gateStatus.value = `Dart ${dartGate.isValid} · u2 ${gateForm.isValid}`);
    }
    const gateParam = [...gateCall.inputParams.values()][0];
    const gateSub = gateParam.onChanged.subscribe(() => refreshGate());
    scope.own(() => gateSub.unsubscribe());
    const dartGateSub = dartGate.onInputChanged.subscribe(() => refreshGate());
    scope.own(() => dartGateSub.unsubscribe());
    refreshGate();

    const gateName = gateParam.name;
    gateGrid.append(span(''), h3('Dart'), h3('u2'));
    const gateCells = [span(gateName, 'u2demo-code'), dartGate.inputs[0].root,
      gateForm.getInput(gateName)!.root];
    for (const cell of gateCells)
      cell.dataset.row = gateName;
    gateGrid.append(...gateCells);
    gateGrid.append(note(REQUIRED));

    await buildAb(dynCall, dynGrid, toDartDyn, toU2Dyn, tickDyn, new Map<string, string[]>([
      ['city', [DEP_REFRESH, PRUNE]],
      ['model', [PROPAGATE]],
      ['site', [SUGGEST]],
      ['cohort', [DEFAULT_OK]],
      ['seed', [DEFAULT_BROKEN]],
      ['regions', [LIST_DYN]],
    ]));

    await buildAb(w3Call, w3Grid, toDartW3, toU2W3, tickW3, new Map<string, string[]>([
      ['df', [P_ICONS, W3_CLOSED, W3_EXTERNAL]],
      ['num', [W3_ORDER, W3_AUTOPICK, P_COMBO]],
      ['cat', [W3_NOMATCH]],
      ['metrics', [W3_LIST]],
    ]));

    plain.append(Scope.runWith(scope, () => funcForm(demoFunc!.prepare())).root);
    // the Dart form writes computed defaults (and the table auto-fill and column auto-picks)
    // during forFuncCall, before the tick subscriptions exist — re-read the calls once so the
    // `inputs` lines start honest
    tick.value++;
    tickDyn.value++;
    tickW3.value++;
  }

  const page = divV([
    prose(INTRO),
    prose(RUN),
    grid,
    divH([span('sync = '), span(computed(() =>
      `u2 → Dart ${toDart.value} · Dart → u2 ${toU2.value}`))], 'u2demo-status'),
    divH([span('fields = '), span(fields)], 'u2demo-status'),
    divH([span('inputs = '), span(computed(() => {
      tick.value;
      return inputsText(call);
    }))], 'u2demo-status'),
    h3('isValid — required parameter'),
    prose(GATE),
    gateGrid,
    divH([span('isValid = '), span(gateStatus)], 'u2demo-status'),
    h3('Dynamic sources'),
    prose(DYN_INTRO),
    dynGrid,
    divH([span('sync = '), span(computed(() =>
      `u2 → Dart ${toDartDyn.value} · Dart → u2 ${toU2Dyn.value}`))], 'u2demo-status'),
    divH([span('inputs = '), span(computed(() => {
      tickDyn.value;
      return inputsText(dynCall);
    }))], 'u2demo-status'),
    h3('Tables & columns'),
    prose(W3_INTRO),
    w3Grid,
    divH([span('sync = '), span(computed(() =>
      `u2 → Dart ${toDartW3.value} · Dart → u2 ${toU2W3.value}`))], 'u2demo-status'),
    divH([span('inputs = '), span(computed(() => {
      tickW3.value;
      return inputsText(w3Call);
    }))], 'u2demo-status'),
    prose(W3_UNASSOC),
    h3('Deliberate divergences'),
    ...DIVERGENCES.map(prose),
    h3('u2 form, unmodified'),
    prose('How `funcForm` renders without the A/B grid: the same function over its own freshly ' +
      'prepared call — not the shared one above, so the counters stay put — with the shipped ' +
      '`u2-form-category` layout as it ships.'),
    plain,
  ], 'u2demo-page u2demo-wide');

  void build();
  return page;
}

/** A full-width row in the A/B grid carrying a divergence note for the row above it. */
function note(text: string): HTMLElement {
  const el = prose(text);
  el.classList.add('u2demo-ab-note');
  return el;
}
