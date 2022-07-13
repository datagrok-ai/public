import {EXCEL_BLOB_TYPE, FileInput} from './file-input';

export namespace UiUtils {
  export function fileInput(
    initialText = 'Drag-n-drop here',
    initialValue: File | null = null,
    onValueChanged: Function | null = null,
    fileType = EXCEL_BLOB_TYPE
  ) {
    return new FileInput(initialText, initialValue, onValueChanged, fileType);
  }
}
