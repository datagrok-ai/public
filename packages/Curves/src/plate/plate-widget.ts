import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {TYPE} from "datagrok-api/dg";
import {div} from "datagrok-api/ui";

const colorScheme = [DG.Color.white, DG.Color.gray];

export class PlateWidget extends DG.Widget {

  _plateData: DG.DataFrame = DG.DataFrame.create();
  _colorColumn?: DG.Column<number>;
  grid: DG.Grid = DG.Viewer.heatMap(this._plateData);
  _posToRow: Map<String, number> = new Map();
  rows = 0;
  cols = 0;
  colorColumnName: string = '';

  constructor() {
    super(ui.div([], 'curves-plate-widget'));
    this.root.appendChild(this.grid.root);
    this.grid.onCellRender.subscribe((args) => this.renderCell(args));
  }

  static detailedView(table: DG.DataFrame): PlateWidget {
    const pw = new PlateWidget();
    pw.grid.root.style.width = '100%';
    const colorSelector = ui.input.column('Color by', {
      table: table,
      value: pw._colorColumn,
      onValueChanged: (v) => { pw._colorColumn = v; pw.grid.invalidate(); }
    });

    const detailsDiv = div();
    pw.grid.onCurrentCellChanged.subscribe((gc) => {
      ui.empty(detailsDiv);
      const row = pw.dataRow(gc);
      if (row)
        detailsDiv.appendChild(ui.tableFromMap(pw.plateData.rows.get(row).toMap()));
    });

    pw.root.prepend(colorSelector.root);
    pw.root.append(detailsDiv);

    pw.plateData = table;
    return pw;
  }

  get plateData() { return this._plateData; }
  set plateData(t: DG.DataFrame) {
    this._plateData = t;
    this._colorColumn = t.columns.firstWhere((col) => col.type == TYPE.FLOAT);
    let rowCol: DG.Column<number> = this._plateData.col('row')!;
    let colCol: DG.Column<number> = this._plateData.col('col')!;
    for (let i = 0; i < this._plateData.rowCount; i++)
      this._posToRow.set(`${rowCol.get(i)}:${colCol.get(i)}`, i);

    this.rows = rowCol.stats.max;
    this.cols = colCol.stats.max;

    // row header + all columns
    this.grid.dataFrame = DG.DataFrame.create(this.rows);
    this.grid.columns.clear();
    for (let i = 0; i <= 12; i++)
      this.grid.columns.add({gridColumnName: i.toString(), cellType: 'string'});

    this.grid.invalidate();
  }

  dataRow(gc: DG.GridCell): number | undefined {
    return this._posToRow.get(`${gc.gridRow + 1}:${gc.gridColumn.idx}`);
  }

  renderCell(args:DG.GridCellRenderArgs) {
    const gc = args.cell;
    args.g.fillStyle = 'grey'; //(args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
    args.g.strokeStyle = 'grey';
    let g = args.g;
    let x = args.bounds.x, y = args.bounds.y, w = args.bounds.width, h = args.bounds.height;
    const dataRow  = this.dataRow(gc);
    g.textAlign = 'center';
    g.textBaseline = 'middle';
    g.font = `${Math.ceil(Math.min(...[16, w - 1, h - 1]))}px  Roboto, Roboto Local`;

    // column header
    if (gc.isColHeader && gc.gridColumn.idx > 0)
      g.fillText('' + gc.gridColumn.idx, x + w / 2, y + h / 2);
    // row header
    else if (gc.gridColumn.idx == 0 && gc.gridRow > 0)
      g.fillText(String.fromCharCode(65 + gc.gridRow), x + w / 2, y + h / 2);
    else if (h > 0 && dataRow) {
      g.beginPath();
      const r = Math.min(h / 2, w / 2) * 0.8;
      g.ellipse(x + w / 2, y + h / 2, r, r, 0, 0, 2 * Math.PI);

      if (this._colorColumn) {
        const color = this._colorColumn.isNumerical
          ? DG.Color.scaleColor(this._colorColumn.get(dataRow)!, this._colorColumn.min, this._colorColumn.max, undefined, colorScheme)
          : this._colorColumn.meta.colors.getColor(dataRow);
        g.fillStyle = DG.Color.toHtml(color);
        g.fill();
      }

      g.stroke();
    }
    args.preventDefault();
  }
}