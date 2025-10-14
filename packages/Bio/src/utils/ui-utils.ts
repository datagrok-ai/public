import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export function getMacromoleculeColumns(): DG.Column<string>[] {
  const columns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (columns === null) {
    grok.shell.error('Current table does not contain macromolecules');
    return [];
  }
  return columns;
}

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export function adjustGridcolAfterRender(
  grid: DG.Grid, colName: string, width: number, rowHeight?: number, dontWait?: boolean
) {
  const update = () => {
    const col = grid.col(colName);
    if (col)
      col.width = width;
    if (rowHeight)
      grid.props.rowHeight = rowHeight;
  };
  if (dontWait) {
    update();
    return;
  }

  const sub = grid.onAfterDrawOverlay.subscribe(() => {
    sub.unsubscribe();
    setTimeout(() => {
      update();
    });
  });
}
