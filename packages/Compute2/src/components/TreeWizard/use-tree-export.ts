import * as grok from 'datagrok-api/grok';
import * as Vue from 'vue';
import {isFuncCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {ExportOption, TreeExportContext, exportWorkflow, workflowExports} from '../../run-exports';
import {downloadFile} from './tree-export';
import type {DriverAPI} from '../../composables/use-reactive-tree-driver';

export function useTreeExport({treeState, states, currentMetaCallData, hasNotSavedEdits}: DriverAPI) {
  const context = (): TreeExportContext => ({
    states: {
      callInfoStates: states.calls,
      validationStates: states.validations,
      consistencyStates: states.consistency,
      pipelineValidations: states.pipelineValidations,
      descriptions: states.descriptions,
    },
    metaStates: states.meta,
    meta: currentMetaCallData.value,
    hasNotSavedEdits: hasNotSavedEdits.value,
  });

  const exports = Vue.computed(() =>
    !treeState.value || isFuncCallState(treeState.value) ? [] : workflowExports(treeState.value));

  const exportHandler = async ({id}: ExportOption) => {
    if (!treeState.value)
      return;
    try {
      downloadFile(await exportWorkflow(treeState.value, id, context()));
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
