import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

function getChoices(gridCell: DG.GridCell): string[] | null {
  const choicesStr = gridCell.tableColumn?.getTag(DG.TAGS.CHOICES);
  if (!choicesStr)
    return null;

  return JSON.parse(choicesStr);
}

/**
 * Renders a comma-separated string value as checkboxes with options retrieved
 * from the column's `.choices` tag.
 * */
export class MultiChoiceCellRenderer extends DG.GridCellRenderer {
  get name() { return 'test'; }

  get cellType() { return 'MultiChoice'; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    const checkEmpty = '\uf0c8';
    const checkSquare = '\uf14a';

    const choices = getChoices(gridCell);
    if (!choices)
      return;
    const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());


    for (let i = 0; i < choices.length; i++) {
      const choice = choices[i];
      const checked = !!values.find((x) => x === choice);
      g.font = '100 14px "Font Awesome 5 Pro"';
      g.fillStyle = DG.Color.toHtml(DG.Color.setAlpha(gridCell.grid.props.cellTextColor, checked ? 255: 150));
      g.fillText(checked ? checkSquare : checkEmpty, x + 4, y + 16 + i * 16);

      g.font = cellStyle.font;
      g.fillText(choice, x + 20, y + 16 + i * 16);
    }
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    const idx = Math.floor((e.offsetY - gridCell.bounds.top - 2) / 16);
    const choices = getChoices(gridCell);
    if (!choices || idx < 0 || idx >= choices.length)
      return;

    const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());
    const itemIdx = values.indexOf(choices[idx]);
    if (itemIdx != -1)
      values.splice(itemIdx, 1)
    else
      values.push(choices[idx]);

    gridCell.setValue(values.join(', '), true);
  }
}
