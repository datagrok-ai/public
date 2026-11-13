import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {take} from 'rxjs/operators';
import {BehaviorSubject} from 'rxjs';
import {FittingView, TargetDescription} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {buildDiffGrokFromFunc} from '@datagrok-libraries/compute-utils/function-views/src/fitting/diff-studio/diff-grok-meta';
import {
  RangeDescription, SensitivityAnalysisView,
} from '@datagrok-libraries/compute-utils/function-views/src/sensitivity-analysis-view';
import {pinView as pinViewHelper} from '../../utils';

export function useRfvAnalyses(
  currentCall: Readonly<Vue.Ref<DG.FuncCall>>,
  callMeta: Readonly<Vue.Ref<Record<string, BehaviorSubject<any>> | undefined>>,
  dgView: DG.ViewBase,
  onFitted: (call: DG.FuncCall) => void,
) {
  const isFittingActive = Vue.ref(false);

  const getRanges = (specificRangeName: string) => {
    const ranges: Record<string, RangeDescription> = {};
    const currentMeta = callMeta.value;
    for (const inputParam of currentCall.value.inputParams.values()) {
      const meta$ = currentMeta?.[inputParam.name];
      const range: RangeDescription = {...(meta$?.value?.[specificRangeName] ?? meta$?.value?.['range'] ?? {})};
      if (range.default == null)
        range.default = inputParam.value;
      ranges[inputParam.name] = range ?? {};
    }
    return ranges;
  };

  const getTargets = () => {
    const targets: Record<string, TargetDescription> = {};
    const currentMeta = callMeta.value;
    for (const outputParam of currentCall.value.outputParams.values()) {
      const meta$ = currentMeta?.[outputParam.name];
      const target: RangeDescription = {...(meta$?.value?.['targetFitting'] ?? {})};
      if (target.default == null)
        target.default = outputParam.value;
      targets[outputParam.name] = target ?? {};
    }
    return targets;
  };

  const pinView = () => pinViewHelper(dgView);

  const runSA = async () => {
    try {
      pinView();
      const ranges = getRanges('rangeSA');
      const diffGrok = await buildDiffGrokFromFunc(currentCall.value.func);
      const inputsLookup = diffGrok?.ivp?.inputsLookup ?? undefined;
      SensitivityAnalysisView.fromEmpty(currentCall.value.func,
        {ranges, diffGrok, inputsLookup, disableLookupDefault: true});
    } catch (e: any) {
      grok.shell.error(e);
    }
  };

  const runFitting = async () => {
    if (isFittingActive.value)
      return;
    pinView();
    isFittingActive.value = true;
    try {
      const currentView = grok.shell.v;
      const ranges = getRanges('rangeFitting');
      const targets = getTargets();
      const diffGrok = await buildDiffGrokFromFunc(currentCall.value.func);
      const inputsLookup = diffGrok?.ivp?.inputsLookup ?? undefined;
      const view = await FittingView.fromEmpty(currentCall.value.func,
        {ranges, targets, acceptMode: true, diffGrok, inputsLookup, disableLookupDefault: true});
      const call = await view.acceptedFitting$.pipe(take(1)).toPromise();
      grok.shell.v = currentView;
      if (call)
        onFitted(call);
    } catch (e: any) {
      grok.shell.error(e);
    } finally {
      isFittingActive.value = false;
    }
  };

  return {isFittingActive: Vue.readonly(isFittingActive), runSA, runFitting};
}
