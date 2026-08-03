import * as grok from 'datagrok-api/grok';
import {
  chemFunctionsDialog as _chemFunctionsDialog,
} from '@datagrok-libraries/statistics/src/compute-functions/dialog';
import {IChemFunctionsDialogResult, IComputeDialogResult, IDescriptorTree,
  HitTriageTemplate} from '../types';
import '../../../css/hit-triage.css';
import {HitAppBase} from '../hit-app-base';
import {HTFunctionOrderingLSKey} from '../consts';

export async function chemFunctionsDialog(app: HitAppBase<any>,
  onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  template: Omit<HitTriageTemplate, 'dataSourceType'>, dialog?: boolean,
  // When true, each function/script/query gets a "re-run when campaign opens" checkbox whose state
  // is persisted in the campaign's compute config (see HitDesignApp.reComputeOnCampaignOpen).
  rerunOnOpen?: boolean,
): Promise<IChemFunctionsDialogResult> {
  return _chemFunctionsDialog(
    app.computeFunctions,
    onOk, onCancel,
    template,
    dialog,
    () => grok.chem.descriptorsTree() as Promise<IDescriptorTree>,
    HTFunctionOrderingLSKey,
    undefined,
    rerunOnOpen ? {
      label: 'Re-run when campaign opens',
      tooltip: 'Recompute this function for the whole molecule column every time this campaign is ' +
        'opened. Useful for database lookups that may resolve after the molecule was first drawn.',
    } : undefined,
  );
}
