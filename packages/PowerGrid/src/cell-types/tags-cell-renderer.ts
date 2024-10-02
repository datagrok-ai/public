// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


@grok.decorators.cellRenderer({
  name: 'Tags',
  cellType: 'Tags',
})
/**
 * Renders a comma-separated string value as checkboxes with options retrieved
 * from the column's `.choices` tag.
 * */
export class TagsCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Tags'; }

  get cellType() { return 'Tags'; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    // somehow dart comes null here, will need to investigate it and fix it, now just a workaround
    if (!gridCell.gridColumn.dart)
      return;
    const values: string[] = gridCell.cell.valueString.split(',').map((s) => s.trim());

    let cx = 2;
    let cy = 3;
    g.textBaseline = 'top';

    const getColor = (tag: string) => {
      const colors = gridCell.gridColumn.temp['catColors'] ??= {};
      if (colors[tag] || (colors[tag] === 0))
        return colors[tag];

      const keys = Object.keys(colors);
      colors[tag] ??= DG.Color.getCategoricalColor(keys.length);
    };

    for (const tag of values) {
      const width = g.measureText(tag).width;
      const drawTag = () => {
        const color = getColor(tag);
        g.fillStyle = DG.Color.toHtml(color);
        g.roundRect(x + cx, y + cy, width + 4, 16, 4);
        g.fill();

        g.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(color));
        g.fillText(tag, x + cx + 2, y + cy + 1);
      };

      if (cx + width <= w) {
        drawTag();
        cx += width + 8;
      }
      else {
        cx = 2;
        cy += 16;
        drawTag();
      }
    }
  }
}
