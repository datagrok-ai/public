// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {GridColumn} from 'datagrok-api/dg';

export function getChoices(column: DG.Column): string[] | null {
  const choices = column?.meta.choices;
  if (!choices || choices.length === 0)
    return null;

  return choices;
}


@grok.decorators.cellRenderer({
  name: 'Multi Choice',
  cellType: 'MultiChoice',
})
/**
 * Renders a comma-separated string value as checkboxes with options retrieved
 * from the column's `.choices` tag.
 * */
export class MultiChoiceCellRenderer extends DG.GridCellRenderer {
  get name() { return 'MultiChoice'; }

  get cellType() { return 'MultiChoice'; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const checkEmpty = '\uf0c8';
    const checkSquare = '\uf14a';

    // somehow dart comes null here, will need to investigate it and fix it, now just a workaround
    if (!gridCell.gridColumn.dart)
      return;
    const choices = getChoices(gridCell.tableColumn!);
    if (!choices)
      return;
    const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());

    for (let i = 0; i < choices.length; i++) {
      const choice = choices[i];
      const checked = !!values.find((x) => x === choice);
      g.font = '100 14px "Font Awesome 5 Pro"';
      g.fillStyle = DG.Color.toHtml(checked ? gridCell.grid.props.cellTextColor : DG.Color.lightGray);
      g.fillText(checked ? checkSquare : checkEmpty, x + 4, y + 16 + i * 16);

      g.font = cellStyle.font;
      g.fillText(choice, x + 20, y + 16 + i * 16);
    }
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.gridColumn.editable)
      return;
    const idx = Math.floor((e.offsetY - gridCell.bounds.top - 2) / 16);
    const choices = getChoices(gridCell.tableColumn!);
    if (!choices || idx < 0 || idx >= choices.length)
      return;

    const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());
    const itemIdx = values.indexOf(choices[idx]);
    if (itemIdx != -1)
      values.splice(itemIdx, 1);
    else
      values.push(choices[idx]);

    gridCell.setValue(values.join(', '), true);
  }

  getDefaultSize(gridColumn: GridColumn): {width?: number | null, height?: number | null} {
    const choices = getChoices(gridColumn.column!);
    if (!choices)
      return {width: 20, height: 20};

    const g = gridColumn.grid.canvas.getContext('2d')!;
    const maxWidth = Math.max(...choices?.map((c) => g.measureText(c).width));
    return {
      width: Math.min(200, maxWidth + 30),
      height: Math.min(200, choices.length * 16 + 4),
    };
  }
}
