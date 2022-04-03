import * as DG from 'datagrok-api/dg';

export function updateDivInnerHTML(div: HTMLElement, content: any) {
  div.innerHTML = '';
  div.append(content);
}

export function createFilters(df: DG.DataFrame) {
  return DG.Viewer.fromType('Filters', df, {
    'showContextMenu': false,
  });
}
