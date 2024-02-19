import * as DG from 'datagrok-api/dg';
import {EXCEL_BLOB_TYPE, FileInput} from './src/file-input';
import {HistoryInput, HistoryInputJSON} from './src/history-input';
import {HistoryPanel} from './src/history-panel';

export namespace UiUtils {
  export function fileInput(
    initialText = 'Drag-n-drop here',
    initialValue: File | null = null,
    onValueChanged: Function | null = null,
    fileType: string | null = EXCEL_BLOB_TYPE,
  ) {
    return new FileInput(initialText, initialValue, onValueChanged, fileType);
  }

  export function historyInput(
    label: string,
    funcName: string,
    stringValueFunc: (currentRun: DG.FuncCall) => string,
    visibleColumnsForGrid: Record<string, (currentRun: DG.FuncCall) => string>,
    visibleColumnsForFilter: string[] = [],
    includeParams = true,
  ) {
    return new HistoryInput(
      label, funcName, stringValueFunc, visibleColumnsForGrid, visibleColumnsForFilter, includeParams);
  }

  export function historyInputJSON(
    label: string,
    funcName: string,
    stringValueFunc: (currentRun: DG.FuncCall) => string,
    visibleColumnsForGrid: Record<string, (currentRun: DG.FuncCall) => string>,
    visibleColumnsForFilter: string[] = [],
    includeParams = true,
    skipDfLoad = false,
  ) {
    return new HistoryInputJSON(
      label, funcName, stringValueFunc, visibleColumnsForGrid, visibleColumnsForFilter, includeParams, skipDfLoad);
  }

  export function historyPanel(objFunc: DG.Func) {
    return new HistoryPanel(objFunc);
  }
}
