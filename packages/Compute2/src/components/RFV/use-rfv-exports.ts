import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {canUseResults} from '../../utils';
import {saveCallToProject, DfExportEntry} from '../../project-export';
import {downloadFile} from '../TreeWizard/tree-export';
import {DEFAULT_EXPORT, callExports, exportCall} from '../../run-exports';

interface ExportItem {
  id: string,
  name: string,
  handler: (arg?: any) => void,
}

export function useRfvExports(currentCall: Readonly<Vue.Ref<DG.FuncCall>>, deps: {
  hideDefaultExport: () => boolean,
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

  const runExport = async (exportId: string) => {
    try {
      downloadFile(await exportCall(currentCall.value, exportId, {
        viewersHook: viewersHook.value,
        metaState: callMeta.value,
        validation: validationState.value,
        consistency: consistencyState.value,
        isOutputOutdated: isOutputOutdated.value,
        completed: isOutputOutdated.value === false,
      }));
    } catch (e: any) {
      grok.shell.error(e);
    }
  };

  const exports = Vue.computed<ExportItem[]>(() => [
    ...callExports(currentCall.value.func, deps.hideDefaultExport())
      .map(({id, name}) => ({id, name, handler: () => runExport(id)})),
    {
      id: 'project',
      name: 'Save as project...',
      handler: () => saveCallToProject(currentCall.value, deps.collectDfExportEntries()),
    },
  ]);

  // The Excel report (or a sole export) also gets a direct ribbon icon next to the menu
  const directExport = Vue.computed(() => {
    const list = exports.value;
    return list.find((e) => e.id === DEFAULT_EXPORT) ?? (list.length === 1 ? list[0] : null);
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
