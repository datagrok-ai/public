import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import {
  isFuncCallState,
  PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {CustomExport, ExportCbInput, ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {disposeViewers, getViewers} from '../../utils';
import {applyCustomExport, reportTree, getExportSummary, reportSummary} from './tree-export';
import type {DriverAPI} from '../../composables/use-reactive-tree-driver';

export function useTreeExport({treeState, states, currentMetaCallData, hasNotSavedEdits}: DriverAPI) {
  const exportStates = () => ({
    callInfoStates: states.calls,
    validationStates: states.validations,
    consistencyStates: states.consistency,
    pipelineValidations: states.pipelineValidations,
    descriptions: states.descriptions,
  });

  const exports = Vue.computed(() => {
    if (!treeState.value || isFuncCallState(treeState.value))
      return [];
    const activeExports: CustomExport[] = [];
    if (!treeState.value.disableDefaultExport) {
      const defaultExport: CustomExport = {
        id: 'default',
        friendlyName: 'Default Excel',
        handler: () => reportTree({
          startDownload: true,
          treeState: treeState.value!,
          meta: currentMetaCallData.value,
          ...exportStates(),
          metaStates: states.meta,
          hasNotSavedEdits: hasNotSavedEdits.value,
        }),
      };
      activeExports.push(defaultExport);
    }
    activeExports.push(...(treeState.value.customExports ?? []));
    return activeExports;
  });

  const exportHandler = async (exportData: CustomExport) => {
    if (!treeState.value)
      return;
    const utils = {
      reportFuncCallExcel: async (
        fc: DG.FuncCall, uuid: string, viewersHook?: ViewersHook, metaState?: Record<string, BehaviorSubject<any>>,
      ) => {
        const viewers = await getViewers(fc, viewersHook, metaState);
        return richFunctionViewReport(
          'Excel',
          fc.func,
          fc,
          viewers,
          states.validations?.[uuid],
          states.consistency?.[uuid],
        ).finally(() => disposeViewers(viewers));
      },
      reportStateExcel: async (state: PipelineState, cb?: (input: ExportCbInput) => Promise<void>) => {
        return reportTree({
          startDownload: false,
          treeState: state,
          meta: currentMetaCallData.value,
          ...exportStates(),
          metaStates: states.meta,
          hasNotSavedEdits: hasNotSavedEdits.value,
          cb,
        });
      },
      getExportSummary: (state: PipelineState) => getExportSummary(state, exportStates()),
      reportSummaryExcel: (state: PipelineState) => reportSummary(getExportSummary(state, exportStates())),
      getFuncCallCustomExports: (fc: DG.FuncCall) => {
        return Utils.getCustomExports(fc.func).map((x) => x.name);
      },
      runFuncCallCustomExport: async (fc: DG.FuncCall, uuid: string, exportName: string) => {
        return applyCustomExport(fc, exportName, {
          startDownload: false,
          validationState: states.validations?.[uuid],
          consistencyState: states.consistency?.[uuid],
          isOutputOutdated: states.calls?.[uuid]?.isOutputOutdated,
          runError: states.calls?.[uuid]?.runError,
        });
      },
    };
    try {
      return await exportData.handler(treeState.value!, utils);
    } catch (e: any) {
      grok.shell.error(e);
    }
  };

  const isTreeReportable = Vue.computed(() => {
    return Object.values(states.calls)
      .map((state) => state?.isOutputOutdated)
      .every((isOutdated) => isOutdated === false);
  });

  return {exports, exportHandler, isTreeReportable};
}
