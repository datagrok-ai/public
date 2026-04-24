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
): Promise<IChemFunctionsDialogResult> {
  return _chemFunctionsDialog(
    app.computeFunctions,
    onOk, onCancel,
    template,
    dialog,
    () => grok.chem.descriptorsTree() as Promise<IDescriptorTree>,
    HTFunctionOrderingLSKey,
  );
}
