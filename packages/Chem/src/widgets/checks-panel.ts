/**
 * The Checks section of the context panel: the foot-gun flags of the current cell, a molecule or a reaction, as
 * utils/footguns/rules.ts words them, the worst first. It shows in Crux mode only, and only for a value with something
 * to say (RDKit mode stays as it was): package.ts's condition asks {@link checksShown}.
 */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Flag, LEVEL_ORDER, Level, checkValue} from '../utils/footguns/checks';

/** Each level's colour and icon (Font Awesome), as the dialog footer's slot shows them. */
const LEVEL_COLOUR: {[level in Level]: string} = {block: 'var(--red-3)', warn: 'var(--orange-2)', note: 'var(--grey-5)'};
const LEVEL_ICON: {[level in Level]: string} = {block: 'ban', warn: 'exclamation-triangle', note: 'info-circle'};
/** The list's class: what the tests and the BDD bindings find it by; each line carries its level and rule. */
export const CHECKS_LIST_CLASS = 'chem-footgun-checks';

/** Whether the Checks section shows for `value` of `semType` (the panel's condition sees the value only; the current
 * object names its cell). */
export function checksShown(value: string, semType: string): boolean {
  return checkValue(value, semType).length > 0;
}

/** The Checks section of a semantic value: one line per flag, the worst first, its level's icon and colour, its words
 * whole. */
export function checksPanel(semValue: DG.SemanticValue): DG.Widget {
  const flags: Flag[] = checkValue(semValue.value, semValue.semType, semValue.cell)
    .sort((a, b) => LEVEL_ORDER[a.level] - LEVEL_ORDER[b.level]);
  const lines = flags.map((f) => {
    const icon = ui.iconFA(LEVEL_ICON[f.level]);
    icon.style.marginRight = '6px';
    const line = ui.divH([icon, ui.divText(f.text)], {style: {color: LEVEL_COLOUR[f.level], alignItems: 'baseline',
      marginBottom: '4px'}});
    line.dataset.level = f.level;
    line.dataset.rule = f.rule;
    return line;
  });
  return new DG.Widget(ui.divV(lines, CHECKS_LIST_CLASS));
}
