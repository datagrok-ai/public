import * as DG from 'datagrok-api/dg';
import {Scope, divV, signal, span} from '@datagrok-libraries/u2';
import {funcForm, functionInput} from '@datagrok-libraries/u2/dg';
import {prose} from './func-convergence';
import {readout} from './pages/common';

const NOT_RUN = '(not run)';
const RUNNING = 'running…';

const INTRO = 'Pick a function on top (the `FunctionInput` control) and its `funcForm` builds ' +
  'below with the two standard history properties on: Run (`showRun`) executes the call and ' +
  'saves it on the server, so it joins the run history with its own run number; the history ' +
  'icon (`showHistory`) opens the saved runs of this function — rows rendered by the ' +
  'platform\'s FuncCall handler — and clicking one copies that run\'s input values back into ' +
  'the form. Run, tweak, run again, then pull any earlier run back from the popup.';

export function funcHistoryPage(): HTMLElement {
  const scope = Scope.ambient!;
  const name = signal('');
  const input = scope.build(() => functionInput('Function', {bind: name, scalarOnly: true}));
  const holder = divV([], 'u2demo-func-history-form');
  const result = signal(NOT_RUN);
  let formScope = new Scope();
  scope.own(() => formScope.dispose());

  scope.effect(() => {
    const n = name.value;
    formScope.dispose();
    formScope = new Scope();
    holder.textContent = '';
    result.value = NOT_RUN;
    if (n === '')
      return;
    const f = DG.Func.find({}).find((x) => x.nqName === n);
    if (f === undefined)
      return;
    try {
      const call = f.prepare();
      const form = formScope.build(() => funcForm(call, {showHistory: true, showRun: true,
        onRun: () => result.value = f.outputs.map((p) => `${p.name}: ${call.outputs[p.name]}`).join(' · ') || 'done'}));
      // the form only greys Run out while the call runs and saves; say it where the result lands
      const busy = new MutationObserver(() => {
        if (form.root.getAttribute('aria-busy') === 'true')
          result.value = RUNNING;
        else if (result.peek() === RUNNING)
          result.value = 'failed — see the error message';
      });
      busy.observe(form.root, {attributes: true, attributeFilter: ['aria-busy']});
      formScope.own(() => busy.disconnect());
      holder.append(form.root);
    }
    catch (e: any) {
      holder.append(span(String(e?.message ?? e), 'u2demo-error'));
    }
  });

  return divV([
    prose(INTRO),
    input.root,
    holder,
    readout('output', result),
  ], 'u2demo-page u2demo-func-history');
}
