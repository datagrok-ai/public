import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {signal, computed, Scope, div, divV, divH, span, h3, button} from '@datagrok-libraries/u2';
import {funcForm, FuncCallForm} from '@datagrok-libraries/u2/src/dg/index.js';
import {ensureFuncs, fmt, prose} from './func-convergence';

const INTRO = 'One FuncCall, two editors: the platform\'s `DG.InputForm.forFuncCall` on the ' +
  'left, u2\'s `funcForm` on the right — both edit the same freshly prepared call, so an edit ' +
  'on either side shows up on the other. The showcase function opens preselected: static and ' +
  'dynamic choices, a dependent country/city pair, a DataFrame choice propagating `mpg`/`cyl`, ' +
  'typed-text suggestions, a computed default, bounded numbers, sections. Type in the list on ' +
  'the left (`#tag` terms work) to switch to any platform function; `Run` executes the shared ' +
  'call with whatever the forms hold.';

let showcaseFunc: DG.Func | null = null;

/** Registered once per session (the ensureFuncs pattern): one function whose parameters cover
 * the feature matrix both form generators support cleanly — the known Dart-side traps
 * (friendlyName captions, Radio, `editor:` hints, broken defaults) stay on the bench. */
function ensureShowcase(): void {
  if (showcaseFunc != null)
    return;
  ensureFuncs();
  showcaseFunc = grok.functions.register({
    signature: 'string fceShowcase(string name, string stage, bool active, int replicates, ' +
      'double doseLevel, string country, string city, string model, int mpg, int cyl, ' +
      'string site, int cohort, list regions, datetime started)',
    run: (...args: any[]) => args.map((x) => String(x ?? '')).join(' '),
  });
  const inputs = new Map(showcaseFunc.inputs.map((p) => [p.name, p] as [string, DG.Property]));
  inputs.get('stage')!.choices = ['Discovery', 'Phase I', 'Phase II'];
  const replicates = inputs.get('replicates')!;
  replicates.min = 1;
  replicates.max = 10;
  replicates.defaultValue = 3;
  const dose = inputs.get('doseLevel')!;
  dose.min = 0;
  dose.max = 1000;
  dose.showSlider = true;
  dose.options['units'] = 'mg';
  dose.category = 'Dosing';
  const country = inputs.get('country')!;
  country.options['choices'] = 'fceW2Countries()';
  country.options['descriptions'] = {FR: 'France', DE: 'Germany', US: 'United States'};
  inputs.get('city')!.options['choices'] = 'fceW2Cities';
  const model = inputs.get('model')!;
  model.options['choices'] = 'fceW2CarsDf()';
  model.options['propagateChoice'] = 'all';
  for (const name of ['country', 'city', 'model', 'mpg', 'cyl'])
    inputs.get(name)!.category = 'Study';
  inputs.get('site')!.options['suggestions'] = 'fceW2Suggest';
  inputs.get('cohort')!.options['default'] = '2 + 2';
  inputs.get('regions')!.options['choices'] = 'fceW2Countries()';
  inputs.get('started')!.category = 'Advanced';
  for (const p of showcaseFunc.inputs)
    p.nullable = true;
}

function signature(f: DG.Func): string {
  const inputs = f.inputs.map((p) => `${p.propertyType} ${p.name}`).join(', ');
  const outputs = f.outputs.map((p) => p.propertyType).join(', ');
  return `${f.nqName}(${inputs})${outputs === '' ? '' : ': ' + outputs}`;
}

export function funcsPage(): HTMLElement {
  const scope = Scope.ambient!;
  ensureShowcase();
  const list = div([], 'u2demo-funcs-list');
  const detail = divV([], 'u2demo-funcs-detail');
  let formScope = new Scope();
  scope.own(() => formScope.dispose());
  let generation = 0;

  async function show(f: DG.Func): Promise<void> {
    const gen = ++generation;
    formScope.dispose();
    formScope = new Scope();
    const shown = formScope;
    detail.textContent = '';
    detail.append(h3(f.friendlyName || f.name));
    if (f.description !== '' && f.description != null)
      detail.append(prose(f.description));
    detail.append(span(signature(f), 'u2demo-code'));
    // arbitrary platform functions land here — a param type funcForm has never seen, a run that
    // throws — and the pane must survive any of them, so this is the one boundary that catches
    try {
      const call = f.prepare();
      const dartColumn = divV([h3('Platform form')], 'u2demo-funcs-col');
      const u2Column = divV([h3('u2 funcForm')], 'u2demo-funcs-col');
      detail.append(divH([dartColumn, u2Column], 'u2demo-funcs-ab'));

      // per-column isolation: the Dart form is known to crash on some functions (friendlyName
      // captions, `editor:` hints) — one side failing still leaves the other side working
      let form: FuncCallForm | undefined;
      try {
        form = Scope.runWith(shown, () => funcForm(call));
        u2Column.append(form.root);
        if (form.unsupported.length > 0)
          u2Column.append(prose('Parameters not supported by this form yet: ' +
            form.unsupported.map((name) => `\`${name}\``).join(', ')));
      }
      catch (e: any) {
        u2Column.append(span(String(e?.message ?? e), 'u2demo-error'));
      }

      Scope.runWith(shown, () => {
        const result = signal('');
        const failed = signal(false);
        const run = button('Run', async () => {
          try {
            await call.call();
            failed.value = false;
            const v = call.getOutputParamValue();
            result.value = v == null ? '(no output)' : fmt(v);
          }
          catch (e: any) {
            failed.value = true;
            result.value = 'Run failed: ' + String(e?.message ?? e);
          }
        }, {primary: true});
        if (form === undefined) {
          run.disabled = true;
          run.title = 'The u2 form failed to build';
        }
        else if (form.inputs.length === 0 && form.unsupported.length > 0) {
          run.disabled = true;
          run.title = 'Nothing to run with: this form can\'t edit ' +
            form.unsupported.map((name) => `'${name}'`).join(', ');
        }
        else {
          const gate = form;
          const required = [...call.inputParams.values()].filter((p) =>
            p.property.isOptional !== true && gate.getInput(p.name) !== undefined);
          const blocked = computed(() => gate.validity.value !== null || required.some((p) => {
            const v = gate.getInput(p.name)!.value.value;
            return Array.isArray(v) ? v.length === 0 : (v == null || v === '');
          }));
          Scope.ambient!.effect(() => {
            run.disabled = blocked.value;
            run.title = blocked.value ? 'Fix the highlighted fields to run' : '';
          });
        }
        const status = divH([span(computed(() => failed.value ? '' : 'result = ')), span(result)],
          'u2demo-status');
        Scope.ambient!.effect(() => status.classList.toggle('u2demo-error', failed.value));
        detail.append(divH([run], 'u2demo-row'), status);
      });

      try {
        const dartForm = await DG.InputForm.forFuncCall(call, {twoWayBinding: true});
        if (gen !== generation)
          return;
        dartColumn.append(dartForm.root);
      }
      catch (e: any) {
        if (gen !== generation)
          return;
        dartColumn.append(span(String(e?.message ?? e), 'u2demo-error'));
      }
    }
    catch (e: any) {
      detail.append(span(String(e?.message ?? e), 'u2demo-hint'));
    }
  }

  void show(showcaseFunc!);

  void (async () => {
    const widget = await DG.Func.byName('FunctionsWidget').apply({}) as DG.FunctionsWidget;
    list.append(widget.root);
    const sub = widget.onActionClicked.subscribe((f: DG.Func) => void show(f));
    scope.own(() => sub.unsubscribe());
    scope.own(() => widget.detach());
  })().catch((e: any) => list.append(span(String(e?.message ?? e), 'u2demo-hint')));

  return divV([
    prose(INTRO),
    divH([list, detail], 'u2demo-funcs-split'),
  ], 'u2demo-page u2demo-wide u2demo-funcs');
}
