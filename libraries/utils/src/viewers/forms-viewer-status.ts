import * as DG from 'datagrok-api/dg';
import type {FormsViewer} from './forms-viewer';

type Box = {x: number, y: number, width: number, height: number};

const CARD = '.d4-multi-form-form:not(.temp)';

/** The cards of the virtual list in the order the user reads them: it lays them out left to right,
 * and `refreshItem` re-appends the card it rebuilt, so the DOM order is not the position order. */
function listCards(v: FormsViewer): HTMLElement[] {
  const cards = Array.from(v.virtualView.root.querySelectorAll(CARD)) as HTMLElement[];
  return cards.sort((a, b) => {
    const ra = a.getBoundingClientRect();
    const rb = b.getBoundingClientRect();
    return (ra.left - rb.left) || (ra.top - rb.top);
  });
}

/** The table row each card position shows: the current row, the mouse-over row, then the records. */
function cardRows(v: FormsViewer): number[] {
  const rows: number[] = [];
  if (v.showCurrentRow)
    rows.push(v.dataFrame.currentRowIdx);
  if (v.showMouseOverRow)
    rows.push(v.dataFrame.mouseOverRowIdx);
  for (const i of Array.from(v.indexes))
    rows.push(i);
  return rows;
}

function fieldsOf(card: HTMLElement): HTMLElement[] {
  return Array.from(card.querySelectorAll('[column]')) as HTMLElement[];
}

function fieldText(el: HTMLElement): string {
  if (el instanceof HTMLCanvasElement)
    return 'canvas';
  const input = el as HTMLInputElement;
  if (input.type === 'checkbox')
    return String(input.checked);
  return input.value ?? '';
}

function backgroundOf(el: HTMLElement): string {
  const color = getComputedStyle(el).backgroundColor;
  const rgb = color.match(/^rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*([\d.]+))?\)$/);
  if (!rgb || (rgb[4] !== undefined && Number(rgb[4]) === 0))
    return '';
  const hex = (s: string) => Number(s).toString(16).padStart(2, '0');
  return `#${hex(rgb[1])}${hex(rgb[2])}${hex(rgb[3])}`.toUpperCase();
}

/** What the Forms viewer shows, for automation: its panes as parts, every card, field, header label,
 * remove icon and sort indicator as a hit area in CSS px of the viewer's root, and the readings a
 * test compares — nothing here is kept between calls, and nothing outside this file computes it. */
export function formsViewerStatus(v: FormsViewer): DG.IWidgetStatus {
  const origin = v.root.getBoundingClientRect();
  const hitAreas: {[name: string]: Box} = {};
  const values: {[name: string]: number | string | boolean} = {};

  const put = (name: string, el: Element | null | undefined) => {
    if (!el)
      return;
    const r = el.getBoundingClientRect();
    if (r.width > 0 && r.height > 0)
      hitAreas[name] = {x: r.left - origin.left, y: r.top - origin.top, width: r.width, height: r.height};
  };

  const readCard = (card: HTMLElement, label: string) => {
    for (const field of fieldsOf(card)) {
      const name = field.getAttribute('column')!;
      const r = field.getBoundingClientRect();
      put(`field ${name} of ${label}`, field);
      values[`${name} of ${label}`] = fieldText(field);
      values[`width of ${name} of ${label}`] = Math.round(r.width);
      values[`height of ${name} of ${label}`] = Math.round(r.height);
      values[`background of ${name} of ${label}`] = backgroundOf(field);
      values[`align of ${name} of ${label}`] = getComputedStyle(field).textAlign;
      values[`font of ${name} of ${label}`] = field.style.font;
      if (values[`field kind of ${name}`] === undefined)
        values[`field kind of ${name}`] = field instanceof HTMLCanvasElement ? 'canvas' : 'input';
    }
  };

  const cards = listCards(v);
  const rows = cardRows(v);
  let shown = 0;
  for (let i = 0; i < cards.length; i++) {
    const label = `card ${i + 1}`;
    put(label, cards[i]);
    readCard(cards[i], label);
    values[`record of ${label}`] = (rows[i] ?? -1) >= 0 ? rows[i] + 1 : '';
    values[`card kind of ${label}`] = v.showCurrentRow && i === v.currentRowPos ? 'current' :
      v.showMouseOverRow && i === v.mouseOverPos ? 'mouse-over' : 'record';
    if (v.showCurrentRow && i === v.currentRowPos) {
      put('current card', cards[i]);
      readCard(cards[i], 'current card');
    }
    if (v.showMouseOverRow && i === v.mouseOverPos) {
      put('mouse-over card', cards[i]);
      readCard(cards[i], 'mouse-over card');
    }
    if ((rows[i] ?? -1) >= 0)
      shown++;
  }

  const pinned = Array.from(v.pinnedFormsDiv.querySelectorAll(CARD)) as HTMLElement[];
  for (let i = 0; i < pinned.length; i++) {
    const label = `pinned card ${i + 1}`;
    put(label, pinned[i]);
    readCard(pinned[i], label);
    // `renderPinnedForms` appends one card per entry of `pinnedRowIndexes`, in that order
    const row = v.pinnedRowIndexes[i] ?? -1;
    values[`record of ${label}`] = row >= 0 ? row + 1 : '';
    values[`card kind of ${label}`] = 'pinned';
  }

  const headerLabels: string[] = [];
  for (const container of Array.from(v.columnHeadersDiv.children) as HTMLElement[]) {
    const label = container.firstElementChild as HTMLElement | null;
    const name = label?.textContent ?? '';
    if (!name)
      continue;
    headerLabels.push(name);
    put(`label ${name}`, label);
    put(`remove ${name}`, container.querySelector('.grok-icon'));
    put(`sort indicator ${name}`, container.querySelector('.d4-multi-form-column-sort-indicator'));
  }

  const sortColumns = v.getSortByColumns();
  values['cards'] = cards.length;
  values['records shown'] = shown + pinned.length;
  values['pinned records'] = pinned.length;
  values['fields shown'] = v.fieldsColumnNames.length;
  values['fields'] = v.fieldsColumnNames.join(', ');
  values['header labels'] = headerLabels.join(', ');
  values['pinned pane shown'] = v.pinnedFormsDiv.style.display !== 'none';
  values['pinned values'] = v.pinnedRowValues.join(', ');
  values['pinned by'] = v.pinnedRowColumnNames.join(', ');
  values['sort column'] = sortColumns.length > 0 ? sortColumns[0] : '';
  values['sort direction'] = sortColumns.length === 0 ? '' : (v.getSortByTypes()[0] ? '↑' : '↓');
  values['current record'] = v.dataFrame.currentRowIdx >= 0 ? v.dataFrame.currentRowIdx + 1 : '';
  values['mouse-over record'] = v.dataFrame.mouseOverRowIdx >= 0 ? v.dataFrame.mouseOverRowIdx + 1 : '';

  return {
    parts: {root: v.root, header: v.columnHeadersDiv, list: v.virtualView.root, pinned: v.pinnedFormsDiv},
    hitAreas, values, shortcuts: {}, events: [], description: null, error: null,
  };
}
