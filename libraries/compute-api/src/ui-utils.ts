/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {UiUtils} from '@datagrok-libraries/compute-utils/shared-components';

//@ts-ignore
export const fileInput = (window.compute.fileInput) as
  ((...args: Parameters<typeof UiUtils.fileInput>) => ReturnType<typeof UiUtils.fileInput>);

//@ts-ignore
export const historyInput = (window.compute.historyInput) as
  ((...args: Parameters<typeof UiUtils.historyInput>) => ReturnType<typeof UiUtils.historyInput>);

//@ts-ignore
export const historyInputJSON = (window.compute.historyInputJSON) as
  ((...args: Parameters<typeof UiUtils.historyInputJSON>) => ReturnType<typeof UiUtils.historyInputJSON>);

//@ts-ignore
export const historyPanel = (window.compute.historyPanel) as
  ((...args: Parameters<typeof UiUtils.historyPanel>) => ReturnType<typeof UiUtils.historyPanel>);
