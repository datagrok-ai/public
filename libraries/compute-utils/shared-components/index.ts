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
    // Label placed before next to the input
    label: string,
    // Function to look history for
    funcName: string,
    options?: {
      // FuncCall props (inputs, outputs, options) to be visible. By default, all input params will be visible.
      // You may add outputs and/or options.
      mainProps?: string[],
      // Custom mapping between prop name and it's extraction logic
      propFuncs?: Record<string, (currentRun: DG.FuncCall) => string>,
      // Custom logic for input's stringValue. By default, mainParams will be used to generate string value
      stringValueFunc?: (currentRun: DG.FuncCall) => string,
    },
  ) {
    return new HistoryInput(
      label, funcName, options,
    );
  }

  export function historyInputJSON(
    // Label placed before next to the input
    label: string,
    // Function to look history for
    funcName: string,
    options?: {
      // FuncCall props (inputs, outputs, options) to be visible. By default, all input params will be visible.
      // You may add outputs and/or options.
      mainProps?: string[],
      // Custom mapping between prop name and it's extraction logic
      propFuncs?: Record<string, (currentRun: DG.FuncCall) => string>,
      // Custom logic for input's stringValue. By default, mainParams will be used to generate string value
      stringValueFunc?: (currentRun: DG.FuncCall) => string,
    },
  ) {
    return new HistoryInputJSON(
      label, funcName, options,
    );
  }

  export function historyPanel(objFunc: DG.Func) {
    return new HistoryPanel(objFunc);
  }
}
