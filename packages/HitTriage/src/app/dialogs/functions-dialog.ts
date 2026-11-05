import * as grok from 'datagrok-api/grok';
import {
  chemFunctionsDialog as _chemFunctionsDialog,
} from '@datagrok-libraries/statistics/src/compute-functions/dialog';
import {IChemFunctionsDialogResult, IComputeDialogResult, IDescriptorTree, TemplateCompute} from '../types';
import '../../../css/hit-triage.css';
import {HitAppBase} from '../hit-app-base';
import {HTFunctionOrderingLSKey} from '../consts';

/** The checkbox that flags a function to be calculated again every time a campaign opens. */
export const rerunOnOpenOption = {
  label: 'Re-run when campaign opens',
  tooltip: 'Recompute this function for the whole molecule column every time this campaign is ' +
    'opened. Useful for database lookups that may resolve after the molecule was first drawn.',
};

/** The descriptors a compute editor offers: the platform's chemistry descriptor groups. */
export const descriptorsTree = () => grok.chem.descriptorsTree() as Promise<IDescriptorTree>;

export async function chemFunctionsDialog(app: HitAppBase<any>,
  onOk: (result: IComputeDialogResult) => void, onCancel: () => void,
  template: {compute: TemplateCompute}, dialog?: boolean,
  // When true, each function/script/query gets a "re-run when campaign opens" checkbox whose state
  // is persisted in the campaign's compute config (see HitDesignApp.reComputeOnCampaignOpen).
  rerunOnOpen?: boolean,
): Promise<IChemFunctionsDialogResult> {
  return _chemFunctionsDialog(
    app.computeFunctions,
    onOk, onCancel,
    template,
    dialog,
    descriptorsTree,
    HTFunctionOrderingLSKey,
    undefined,
    rerunOnOpen ? rerunOnOpenOption : undefined,
  );
}
