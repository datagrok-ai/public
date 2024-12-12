import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import '../../css/chem.css';

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export function addCopyIcon(element: Object | HTMLElement, paneName: string) {
  const copyIcon = ui.icons.copy(() => {}, 'Copy');
  copyIcon.onclick = (e) => {
    let text: string;
    if (element instanceof HTMLElement)
      text = element.innerText;
    else {
      let tableString = '';
      for (const [key, value] of Object.entries(element))
        tableString += `${key}\t${(value as any).innerText}\n`;
      text = tableString;
    }
    navigator.clipboard.writeText(text);
    grok.shell.info('Copied to clipboard');
    e.stopImmediatePropagation();
  };
  copyIcon.classList.add('copy-icon');
  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === paneName) {
      const pane = accPanes[i];
      pane.append(copyIcon);
      pane.parentElement?.addEventListener('mouseenter', () => {copyIcon.style.visibility = 'visible';});
      pane.parentElement?.addEventListener('mouseleave', () => {copyIcon.style.visibility = 'hidden';});
    }
  }
}


export function pickTextColorBasedOnBgColor(bgColor: string, lightColor: string, darkColor: string) {
  const color = (bgColor.charAt(0) === '#') ? bgColor.substring(1, 7) : bgColor;
  const r = parseInt(color.substring(0, 2), 16); // hexToR
  const g = parseInt(color.substring(2, 4), 16); // hexToG
  const b = parseInt(color.substring(4, 6), 16); // hexToB
  return (((r * 0.299) + (g * 0.587) + (b * 0.114)) > 186) ?
    darkColor : lightColor;
}

export function getGridCellColTemp<TValue, TTemp>(
  gridCell: DG.GridCell
): [DG.GridColumn | null, DG.Column<TValue>, TTemp] {
  let temp: TTemp | null = null;
  let gridCol: DG.GridColumn | null = null;
  try {
    gridCol = gridCell.dart ? gridCell.gridColumn : null;
    temp = gridCol ? gridCol.temp as TTemp : null;
  } catch { [gridCol, temp] = [null, null]; }

  const tableCol: DG.Column<TValue> = gridCell.cell.column;
  temp = temp ?? tableCol.temp;

  if (!temp)
    throw new Error(`Grid cell renderer back store (GridColumn or Column) not found.`);
  return [gridCol, tableCol, temp];
}

export function getZoomCoordinates(W0: number, H0: number, x1: number, y1: number, x2: number, y2: number) {
  const W1 = Math.abs(x1 - x2);
  const H1 = Math.abs(y1 - y2);
  const scaleW = W0 / W1;
  const scaleH = H0 / H1;
  const scale = Math.min(scaleW, scaleH);
  const W2 = (W0 / scale) * 5;
  const H2 = (H0 / scale) * 5;
  const left = x1 < x2 ? x1 : x2;
  const top = y1 > y2 ? y1 : y2;
  const zoomLeft = (left + W1 / 2) - W2 / 2;
  const zoomRight = zoomLeft + W2;
  const zoomTop = (top - H1 / 2) + H2 / 2;
  const zoomBottom = zoomTop - H2;
  return {zoomLeft: zoomLeft, zoomRight: zoomRight, zoomTop: zoomTop, zoomBottom: zoomBottom};
}

export function resizeGridColsSize(grid: DG.Grid, colsToResize: string[], width: number, height: number) {
  colsToResize.forEach((colName: string) => grid.columns.byName(colName)!.width = width);
  grid.props.rowHeight = height;
}
