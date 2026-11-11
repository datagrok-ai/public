import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {canUseResults, disposeViewers, getViewers} from '../../utils';
import {saveCallToProject, DfExportEntry} from '../../project-export';
import {applyCustomExport, downloadFile} from '../TreeWizard/tree-export';

interface ExportItem {
  name: string,
  handler: (arg?: any) => void,
}

export function useRfvExports(currentCall: Readonly<Vue.Ref<DG.FuncCall>>, deps: {
  isReportEnabled: Readonly<Vue.Ref<boolean>>,
  hideDefaultExport: () => boolean,
  customExports: Readonly<Vue.Ref<{name: string, function: string}[]>>,
  autorun: Readonly<Vue.Ref<boolean>>,
  callState: () => FuncCallStateInfo | undefined,
  isOutputOutdated: Readonly<Vue.Ref<boolean | undefined>>,
  viewersHook: Readonly<Vue.Ref<ViewersHook | undefined>>,
  callMeta: Readonly<Vue.Ref<Record<string, BehaviorSubject<any>> | undefined>>,
  validationState: Readonly<Vue.Ref<Record<string, ValidationResult> | undefined>>,
  consistencyState: Readonly<Vue.Ref<Record<string, ConsistencyInfo> | undefined>>,
  collectDfExportEntries: () => DfExportEntry[],
}) {
  const {isOutputOutdated, viewersHook, callMeta, validationState, consistencyState} = deps;

  const reportHandler = async (exportName: string) => {
    try {
      downloadFile(await applyCustomExport(currentCall.value, exportName, {
        validationState: validationState.value,
        consistencyState: consistencyState.value,
        isOutputOutdated: isOutputOutdated.value,
      }));
    } catch (e: any) {
      grok.shell.error(e);
    }
  };

  const exports = Vue.computed(() => {
    const activeExports: ExportItem[] = [];
    if (deps.isReportEnabled.value && !deps.hideDefaultExport()) {
      const name = 'Default Excel';
      const handler = async () => {
        try {
          const viewers = await getViewers(currentCall.value, viewersHook.value, callMeta.value);
          const [blob] = await richFunctionViewReport(
            'Excel',
            currentCall.value.func,
            currentCall.value,
            viewers,
            validationState.value,
            consistencyState.value,
          ).finally(() => disposeViewers(viewers));
          DG.Utils.download(`${currentCall.value.func.nqName} - ${Utils.getStartedOrNull(currentCall.value) ?? 'Not completed'}.xlsx`, blob);
        } catch (e: any) {
          grok.shell.error(e);
        }
      };
      activeExports.push({name, handler});
    }
    activeExports.push(...deps.customExports.value.filter((x) => x.function && x.name)
      .map((x) => ({...x, handler: () => reportHandler(x.name)})));
    activeExports.push({
      name: 'Save as project...',
      handler: () => saveCallToProject(currentCall.value, deps.collectDfExportEntries()),
    });
    return activeExports;
  });

  // The Excel report (or a sole export) also gets a direct ribbon icon next to the menu
  const directExport = Vue.computed(() => {
    const list = exports.value;
    return list.find((e) => e.name === 'Default Excel') ?? (list.length === 1 ? list[0] : null);
  });

  // On autorun the export icons stay put regardless of isOutputOutdated
  // (no ribbon rebuild/flicker); guard the click so a stale/in-flight run surfaces a shell
  // message instead of exporting partial output.
  const guardedExport = (handler: () => void) => {
    if (canUseResults(deps.callState(), 'exporting'))
      handler();
  };
  const exportsVisible = Vue.computed(() => deps.autorun.value || !isOutputOutdated.value);

  return {exports, directExport, guardedExport, exportsVisible};
}
