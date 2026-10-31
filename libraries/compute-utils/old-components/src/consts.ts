import * as DG from 'datagrok-api/dg';

export const EXPERIMENTAL_TAG = 'experimental';

export const storageName = `ModelStorage`;

export const ID_COLUMN_NAME = 'ID';
export const EXP_COLUMN_NAME = 'Source';
export const FAVORITE_COLUMN_NAME = 'Is favorite';
export const COMPLETE_COLUMN_NAME = 'Is complete';
export const ACTIONS_COLUMN_NAME = 'Delete';
export const AUTHOR_COLUMN_NAME = 'Author';
export const STARTED_COLUMN_NAME = 'Date';
export const TITLE_COLUMN_NAME = 'Title';
export const DESC_COLUMN_NAME = 'Desc.';
export const TAGS_COLUMN_NAME = 'Tags';
export const VERSION_COLUMN_NAME = 'Version';

export type HistoryOptions = {
  // FuncCall props (inputs, outputs, options) to be visible. By default, all input params will be visible.
  // You may add outputs and/or options.
  visibleProps?: string[],
  // Text to show if no runs has been found
  fallbackText?: string,
  // Flag to show per-element edit and delete actions. Default is false
  showActions?: boolean,
  // Flag to show delete and compare actions. Default is false
  showBatchActions?: boolean,
  // Flag used in HistoryPanel. Default is false
  isHistory?: boolean,
  // Custom mapping between prop name and it's extraction logic. Default is empty
  propFuncs?: Record<string, (currentRun: DG.FuncCall) => string>,
  // Versions handing, default none
  version?: string,
  allowOtherVersions?: boolean,
}

export const HISTORY_SUPPORTED_COL_TYPES =
  Object.values(DG.COLUMN_TYPE).filter((type: any) => type !== DG.TYPE.DATA_FRAME);
