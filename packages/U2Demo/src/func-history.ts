import * as DG from 'datagrok-api/dg';
import {Scope, divV, signal, span} from '@datagrok-libraries/u2';
import {funcForm, functionInput} from '@datagrok-libraries/u2/src/dg/index.js';
import {prose} from './func-convergence';

const INTRO = 'Pick a function on top (the `FunctionInput` control) and its `funcForm` builds ' +
  'below with the two standard history properties on: **Run** (`showRun`) executes the call and ' +
  'saves it on the server, so it joins the run history with its own run number; the **history ' +
  'icon** (`showHistory`) opens the saved runs of this function — rows rendered by the ' +
  'platform\'s FuncCall handler — and clicking one copies that run\'s input values back into ' +
  'the form. Run, tweak, run again, then pull any earlier run back from the popup.';

export function funcHistoryPage(): HTMLElement {
  const scope = Scope.ambient!;
  const name = signal('');
  const input = Scope.runWith(scope, () => functionInput('Function', {bind: name, scalarOnly: true}));
  const holder = divV([], 'u2demo-func-history-form');
  let formScope = new Scope();
  scope.own(() => formScope.dispose());

  scope.effect(() => {
    const n = name.value;
    formScope.dispose();
    formScope = new Scope();
    holder.textContent = '';
    if (n === '')
      return;
    const f = DG.Func.find({}).find((x) => x.nqName === n);
    if (f === undefined)
      return;
    try {
      const form = Scope.runWith(formScope, () =>
        funcForm(f.prepare(), {showHistory: true, showRun: true}));
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
  ], 'u2demo-page u2demo-func-history');
}
