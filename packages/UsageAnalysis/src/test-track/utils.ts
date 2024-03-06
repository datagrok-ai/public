// import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {colors} from '../utils';

export {colors};
export const FILENAME = 'test-cases.csv';
export const PASSED = 'passed';
export const FAILED = 'failed';
export const SKIPPED = 'skipped';
export type Status = typeof PASSED | typeof FAILED | typeof SKIPPED;

const map = {
  [PASSED]: {name: 'check', color: 'var(--green-2)'},
  [FAILED]: {name: 'times', color: 'var(--red-3)'},
  [SKIPPED]: {name: 'forward', color: 'var(--orange-2)'},
};

export async function loadFileAsText(name: string): Promise<string> {
  return await _package.files.readAsText(name);
}

export async function readDataframe(tableName: string): Promise<DG.DataFrame> {
  const file = await loadFileAsText(tableName);
  const df = DG.DataFrame.fromCsv(file);
  df.name = tableName.replace('.csv', '');
  return df;
}

export async function writeDataframe(tableName: string, data: string): Promise<void> {
  await _package.files.writeAsText(tableName, data);
}

export function getIcon(name: string, options?: {style?: string, class?: string[], color?: string, id?: string}): HTMLElement {
  const icon = ui.iconFA(name);
  icon.classList.replace('fal', options?.style ?? 'far');
  if (options?.class)
    icon.classList.add(...options.class);
  if (options?.color)
    icon.style.color = options.color;
  if (options?.id)
    icon.id = options.id;
  return icon;
}

export function getStatusIcon(status: typeof PASSED | typeof FAILED | typeof SKIPPED): HTMLElement {
  const obj = map[status];
  return getIcon(obj.name, {class: ['tt-status-icon'], color: obj.color});
}
