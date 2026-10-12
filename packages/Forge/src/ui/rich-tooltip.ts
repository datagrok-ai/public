import * as ui from 'datagrok-api/ui';
import '../../css/forge.css';

/** About two lines of the tooltip's width. */
const MAX_LIST_LENGTH = 90;

interface RichTooltipContent {
  title: string;
  /** The main lines. */
  main: HTMLElement[];
  /** The minor lines, given the grey `forge-tooltip-minor` class by the caller. */
  minor: HTMLElement[];
}

/** Forge's tooltip of an entity: a bold title, the main lines, then the minor lines in grey; at most about 320 px
 * wide. */
export function richTooltip({title, main, minor}: RichTooltipContent): HTMLElement {
  return ui.divV([ui.divText(title, 'forge-tooltip-title'), ...main, ...minor], 'forge-tooltip');
}

/** [list] (a new `ul` by default) emptied and given one bullet per item of [items]. */
export function bulletList(items: (string | HTMLElement)[][], list = ui.element('ul')): HTMLElement {
  ui.empty(list);
  for (const item of items) {
    const li = ui.element('li');
    li.append(...item);
    list.append(li);
  }
  return list;
}

/** [items] joined by commas, cut after about two lines with `+N more`. */
export function shortList(items: string[]): string {
  let text = '';
  for (let i = 0; i < items.length; i++) {
    const next = i === 0 ? items[i] : `${text}, ${items[i]}`;
    if (i > 0 && next.length > MAX_LIST_LENGTH)
      return `${text} +${items.length - i} more`;
    text = next;
  }
  return text;
}
