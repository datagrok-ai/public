import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {funcs} from '../../package-api';

export const sourceRunLabel = (isReadonly: boolean, title?: string) => {
  const run = title ? `'${title}'` : 'a completed';
  return isReadonly ? `Uses ${run} run` : `Based on ${run} run`;
};

export const sourceRunTooltip = (isReadonly: boolean, title?: string) =>
  `${sourceRunLabel(isReadonly, title)}. Click for options`;

export const openSourceRun = (id: string) => funcs.openWorkflowRun(id).catch((e) => grok.shell.error(e));

export interface SourceRunCompareActions {
  compareWithOriginal: () => void;
  addToComparison: () => void;
}

/** Comparison actions are passed for editable items only: a read-only copy equals its original. */
export const showSourceRunMenu = (
  id: string, isReadonly: boolean, causedBy: MouseEvent, title?: string, compare?: SourceRunCompareActions,
) => {
  const menu = DG.Menu.popup()
    .header(sourceRunLabel(isReadonly, title))
    .item(title ? `Open '${title}'` : 'Open original run', () => openSourceRun(id));
  if (compare) {
    menu.item('Compare with original', compare.compareWithOriginal);
    menu.item('Add to comparison...', compare.addToComparison);
  }
  menu.show({causedBy});
};
