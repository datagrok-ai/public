import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  isFuncCallState,
  isParallelPipelineState,
  isSequentialPipelineState, isStaticPipelineState, PipelineState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {zipSync, Zippable} from 'fflate';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';

export function findTreeNode(uuid: string, state: PipelineState): PipelineState | undefined {
  const notVisitedStates = [state];

  while (notVisitedStates.length > 0) {
    const currentState = notVisitedStates.pop()!;

    if (currentState.uuid === uuid)
      return currentState;

    if (
      isParallelPipelineState(currentState) ||
        isSequentialPipelineState(currentState) ||
        isStaticPipelineState(currentState)
    )
      notVisitedStates.push(...currentState.steps);
  }
};

export function findTreeNodeParrent(uuid: string, state: PipelineState): PipelineState | undefined {
  const notVisitedStates = [state];

  while (notVisitedStates.length > 0) {
    const currentState = notVisitedStates.pop()!;

    if (
      isParallelPipelineState(currentState) ||
        isSequentialPipelineState(currentState) ||
        isStaticPipelineState(currentState)
    ) {
      for (const item of currentState.steps) {
        if (item.uuid == uuid)
          return currentState;
      }
      notVisitedStates.push(...currentState.steps);
    }
  }
};


export async function reportStep(treeState?: PipelineState) {
  if (treeState) {
    const zipConfig = {} as Zippable;

    const reportStep = async (state: PipelineState, previousPath: string = '', idx: number = 1) => {
      if (isFuncCallState(state) && state.funcCall) {
        const funccall = state.funcCall;

        const blob = await Utils.richFunctionViewReport(
          'Excel',
          funccall.func,
          funccall,
          Utils.dfToViewerMapping(funccall),
        );

        const validatedFilename = Utils.replaceForWindowsPath(
          `${String(idx).padStart(3, '0')}_${Utils.getFuncCallDefaultFilename(funccall)}`,
        );
        const validatedFilenameWithPath = `${previousPath}/${validatedFilename}`;

        zipConfig[validatedFilenameWithPath] =
                      [new Uint8Array(await blob.arrayBuffer()), {level: 0}];
      }

      if (
        isSequentialPipelineState(state) ||
          isParallelPipelineState(state) ||
          isStaticPipelineState(state)
      ) {
        const nestedPath = `${String(idx).padStart(3, '0')}_${state.friendlyName ?? state.nqName}`;
        let validatedNestedPath = Utils.replaceForWindowsPath(nestedPath);

        if (previousPath.length > 0) validatedNestedPath = `${previousPath}/${validatedNestedPath}`;

        for (const [idx, stepState] of state.steps.entries())
          await reportStep(stepState, validatedNestedPath, idx + 1);
      }
    };

    await reportStep(treeState);

    DG.Utils.download(
      `${treeState.friendlyName ?? treeState.configId}.zip`,
      new Blob([zipSync(zipConfig)]),
    );

  }
}
