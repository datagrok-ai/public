import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {axisBottom, scaleBand, scaleLinear, select, color} from 'd3';
import {ChemPalette} from '../utils/chem-palette';
//@ts-ignore: ignore this shit
import * as rxjs from 'rxjs';
import $ from 'cash-dom';
const cp = new ChemPalette('grok');

export function addViewerToHeader(grid: DG.Grid, viewer: Promise<DG.Widget>) {
  viewer.then((viewer) => {
    const barchart = viewer as StackedBarChart;
    rxjs.fromEvent(grid.overlay, 'mousemove').subscribe((mm:any) => {
      mm = mm as MouseEvent;
      const cell = grid.hitTest(mm.offsetX, mm.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids') {
        barchart.highlight(cell, mm.offsetX, mm.offsetY);
      }
    });

    rxjs.fromEvent(grid.overlay, 'click').subscribe((mm:any) => {
      mm = mm as MouseEvent;

      const cell = grid.hitTest(mm.offsetX, mm.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType == 'aminoAcids') {
        barchart.beginSelection(mm);
        return;
      }
      barchart.unhighlight();
    });
    rxjs.fromEvent(grid.overlay, 'mouseout').subscribe((_: any) => {
      barchart.unhighlight();
    });

    barchart.tableCanvas = grid.canvas;
    grid.setOptions({'colHeaderHeight': 200});
    grid.onCellTooltip((cell, x, y) => {
      if (cell.tableColumn) {
        if (['aminoAcids', 'alignedSequence'].includes(cell.tableColumn.semType) ) {
          if ( !cell.isColHeader) {
            cp.showTooltip(cell, x, y);
            return true;
          } else {
            if (barchart.highlighted) {
              ui.tooltip.show(ui.divV([ui.divText(barchart.highlighted.aaName)]), x, y);
            }
            return true;
          }
        }
      }
    });
    grid.onCellRender.subscribe((args) => {
      args.g.save();
      args.g.beginPath();
      args.g.rect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
      args.g.clip();

      if (args.cell.isColHeader && args.cell.tableColumn?.semType == 'aminoAcids') {
        barchart.renderBarToCanvas(
          args.g,
          args.cell,
          args.bounds.x,
          args.bounds.y,
          args.bounds.width,
          args.bounds.height,
        );
        args.preventDefault();
      }
      args.g.restore();
    });
  });
}


export class StackedBarChart extends DG.JsViewer {
    public dataEmptyAA: string;
    public initialized: boolean;
    highlighted:{'colName':string, 'aaName':string}|null = null;
    private ord: { [Key: string]: number; } = {};
    private margin: { top: number; left: number; bottom: number; right: number } = {
      top: 10,
      right: 10,
      bottom: 50,
      left: 10,
    };
    private yScale: any;
    private xScale: any;
    private data: { 'name': string, 'data': { 'name': string, 'count': number, 'selectedCount': number }[] }[] = [];
    private selectionMode: boolean = false;
    public aminoColumnNames: string[] = [];
    // @ts-ignore

    private aminoColumnIndices: { [Key: string]: number; } = {};
    private aggregatedTables: { [Key: string]: DG.DataFrame; } = {};
    private aggregatedTablesUnselected: { [Key: string]: DG.DataFrame; } = {};
    private max = 0;
    private barStats: { [Key: string]: { 'name': string, 'count': number, 'selectedCount': number }[] } = {};
    tableCanvas: HTMLCanvasElement | undefined;
    private registered: { [Key: string]: DG.GridCell } = {};

    constructor() {
      super();
      this.dataEmptyAA = this.string('dataEmptyAA', '-');
      this.initialized = false;
    }

    init() {
      const groups: [string[], string][] = [
        [['C', 'U'], 'yellow'],
        [['G', 'P'], 'red'],
        [['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'], 'all_green'],
        [['R', 'H', 'K'], 'light_blue'],
        [['D', 'E'], 'dark_blue'],
        [['S', 'T', 'N', 'Q'], 'orange']];

      let i = 0;
      groups.forEach((item) => {
        i++;
        // eslint-disable-next-line guard-for-in
        for (const obj in item[0]) {
          this.ord[item[0][obj]] = i;
        }
      });
      this.yScale = scaleLinear();
      this.xScale = scaleBand();
      this.data = [];

      this.aminoColumnNames = [];
    }

    // Stream subscriptions
    onTableAttached() {
      this.init();
      if (this.dataFrame) {
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
        this.computeData(this.dataFrame);
      }
    }

    // Cancel subscriptions when the viewer is detached
    detach() {
      this.subs.forEach((sub) => sub.unsubscribe());
    }

    computeData(df: DG.DataFrame) {
      this.data = [];
      this.aminoColumnNames = [];
      this.aminoColumnIndices = {};

      df.columns.names().forEach((name: string) => {
        {
          // @ts-ignore
          if (df.getCol(name).semType === 'aminoAcids' &&
                    !df.getCol(name).categories.includes('COOH') &&
                    !df.getCol(name).categories.includes('NH2')) {
            this.aminoColumnIndices[name] = this.aminoColumnNames.length + 1;
            this.aminoColumnNames.push(name);
          }
        }
      });

      this.aggregatedTables = {};
      this.aggregatedTablesUnselected = {};
      const buf1 = df.selection.getBuffer();
      const buf2 = df.filter.getBuffer();
      const resbuf = new Int32Array(df.rowCount);

      for (let i = 0; i < buf2.length; i++) {
        resbuf[i] = buf1[i] & buf2[i];
      }


      const mask = DG.BitSet.fromBytes(resbuf.buffer, df.rowCount);
      if (mask.trueCount !== df.filter.trueCount) {
        this.selectionMode = true;
        this.aminoColumnNames.forEach((name) => {
          this.aggregatedTables[name] = df
            .groupBy([name])
            .whereRowMask(df.filter)
            .add('count', name, `${name}_count`)
            .aggregate();
          const buf1 = df.selection.getBuffer();
          const buf2 = df.filter.getBuffer();
          const resbuf = new Int32Array(df.rowCount);

          for (let i = 0; i < buf2.length; i++) {
            resbuf[i] = buf1[i] & buf2[i];
          }


          // @ts-ignore
          const mask = DG.BitSet.fromBytes(resbuf.buffer, df.rowCount);
          // @ts-ignore
          this.aggregatedTablesUnselected[name] = df
            .groupBy([name])
            .whereRowMask(mask)
            .add('count', name, `${name}_count`)
            .aggregate();
        });
      } else {
        this.selectionMode = false;
        this.aminoColumnNames.forEach((name) => {
          // @ts-ignore
          this.aggregatedTables[name] = df
            .groupBy([name])
            .whereRowMask(df.filter)
            .add('count', name, `${name}_count`)
            .aggregate();
        },
        );
      }
      this.data = [];
      this.barStats = {};
      for (const [name, df] of Object.entries(this.aggregatedTables)) {
        const colObj: {
                'name': string, 'data':
                    { 'name': string, 'count': number, 'selectedCount': number }[]
            } =
                {'name': name, 'data': []};
        this.barStats[colObj['name']] = colObj['data'];
        this.data.push(colObj);
        for (let i = 0; i < df.rowCount; i++) {
          const amino = df.getCol(name).get(i);
          const aminoCount = df.getCol(`${name}_count`).get(i);
          if ((!amino) || amino === this.dataEmptyAA) {
            continue;
          }
          const aminoObj = {'name': amino, 'count': aminoCount, 'selectedCount': 0};
          colObj['data'].push(aminoObj);
          for (let j = 0; j < this.aggregatedTablesUnselected[name].rowCount; j++) {
            const unsAmino = this.aggregatedTablesUnselected[name].getCol(`${name}`).get(j);
            if (unsAmino == amino) {
              aminoObj['selectedCount'] = this.aggregatedTablesUnselected[name]
                .getCol(`${name}_count`)
                .get(j);
              break;
            }
          }
        }
        colObj['data'] = colObj['data'].sort((o1, o2) => {
          if (this.ord[o1['name']] > this.ord[o2['name']]) {
            return -1;
          }
          if (this.ord[o1['name']] < this.ord[o2['name']]) {
            return 1;
          }

          return 0;
        });
      }
      this.max = df.filter.trueCount;
    }

    renderBarToCanvas(g: CanvasRenderingContext2D, cell: DG.GridCell, x: number, y: number, w: number, h: number) {
      const margin = 0.2;
      const innerMargin = 0.02;
      const selectLineRatio = 0.1;
      x = x + w * margin;
      y = y + h * margin / 4;
      w = w - w * margin * 2;
      h = h - h * margin / 2 - w;
      g.fillStyle = 'black';
      g.textBaseline = 'top';
      g.font = `${w / 2}px`;
      // eslint-disable-next-line no-unused-vars

      const name = cell.tableColumn!.name;
      g.fillText(name,
        x+w / 4,
        y + h+ w / 4 );
      const barData = this.barStats[name];
      let sum = 0;
      barData.forEach((obj) => {
        sum += obj['count'];
      });
      let curSum = 0;

      barData.forEach((obj, index) => {
        const sBarHeight = h * obj['count'] / this.max;
        const gapSize = sBarHeight * innerMargin;
        g.fillStyle = cp.getColor(obj['name']);
        g.fillRect(
          x,
          y + h * (this.max - sum + curSum) / this.max + gapSize / 2,
          w,
          sBarHeight - gapSize);
        if (w <= sBarHeight) {
          g.fillStyle = 'rgb(0,0,0)';
          g.font = `${w / 2}px`;
          // eslint-disable-next-line no-unused-vars
          const [_c, aar, _p] = cp.getColorAAPivot(obj['name']);
          g.fillText(aar,
            x + w / 4,
            y + h * (this.max - sum + curSum) / this.max + gapSize / 2 + (sBarHeight - gapSize) / 2 - w / 4);
        }

        if (this.selectionMode && obj['selectedCount'] > 0) {
          g.fillStyle = 'rgb(255,165,0)';
          g.fillRect(
            x - w * selectLineRatio * 1.5,
            y + h * (this.max - sum + curSum) / this.max + gapSize / 2,
            w * selectLineRatio,
            h * obj['selectedCount'] / this.max - gapSize);
        }

        // @ts-ignore
        if (this.dataFrame.currentRow[name] === obj['name']) {
          g.strokeStyle = 'rgb(0,0,0)';
          g.strokeRect(
            x,
            y + h * (this.max - sum + curSum) / this.max + gapSize / 2,
            w,
            sBarHeight - gapSize);
        }

        curSum += obj['count'];
      });
      return;
    }

    render(computeData = true) {
      const df = this.dataFrame!;
      if (computeData) {
        this.computeData(df);
      }
      if (this.tableCanvas) {
        return;
        for (const name of this.aminoColumnNames) {
          this.renderBar(name);
        }
        return;
      }


      // @ts-ignore

      if (!this.root.parentElement) {
        throw new Error('No parent element');
      }
      const
        width = this.root.parentElement.clientWidth;
      this.root.style.width = `${width}px`;
      const
        height = this.root.parentElement.clientHeight;
      this.root.style.width = `${height}px`;
      const
        innerWidth = width - this.margin.left - this.margin.right;
      const
        innerHeight = height - this.margin.top - this.margin.bottom;

      const getColor = cp.getColor;
      const
        scope = this;

      $(this
        .root,
      ).empty();

      const
        svg = select(this.root).append('svg')
          .attr('width', width)
          .attr('height', height);
      svg.attr('style', 'z-index:1');
      const
        g = svg.append('g').attr('transform', `translate(${this.margin.left}, ${this.margin.top})`);
      const
        x = this.xScale
          .domain(this.data.map((d) => d['name']))
          .padding(0.3)
          .align(0.3)
          .rangeRound([0, innerWidth]);

      const
        y = this.yScale
          .domain([0, this.max]).nice()
          .rangeRound([innerHeight, 0]);

      // const aminoScale = scaleBand().domain(Object.entries(this.data[`amino_count`])
      //     .map((d) => (`${d[0]}:${d[1]}`)))
      //     .rangeRound([0, innerHeight]);

      const cAxis = axisBottom(x);
      cAxis.tickSize(0);
      const
        colAxis = g.append('g').call(cAxis)
          .attr('transform', `translate(0, ${innerHeight + 10})`);
      colAxis
        .attr('font-family', '"Open Sans", sans-serif')
        .attr('fill', 'black');
      colAxis.select('.domain').remove();


      // let aAxis = axisLeft(aminoScale)
      // aAxis.tickSize(0);
      //
      // let aminoAxis = g.append("g").call(aAxis)
      // aminoAxis.select(".domain").remove();
      // aminoAxis.attr("font-family", "common sans")
      //     .attr("fill", "black")
      // aminoAxis.selectAll('.tick')
      //     .on('click', (event) => {
      //         scope.dataFrame.selection.handleClick(i => {
      //             let amino_name = event.srcElement.textContent.split(':')[0]
      //             let res = false
      //             this.aminoColumnNames.forEach(name => {
      //                 if (scope.dataFrame.getCol(name).get(i) === amino_name) {
      //                     res = true;
      //                 }
      //             })
      //             return res;
      //         }, event)
      //     });

      g
        .selectAll(
          '.group',
        )
        .data(this

          .data,
        )
        .attr(
          'class'
          ,
          'group',
        )
        .enter()

        .append(
          'g',
        )
        .each(
          function(d: any, i: any) {
            d['data'].map((obj: any, j: any, arr: any) => {
              const xPos = x(d['name']);
              const yPos = ((e) => {
                let sum = 0;
                arr.map((obj1: any, k: any) => {
                  if (k < j) {
                    sum = sum + obj1['count'];
                  }
                });
                return y(obj['count'] + sum);
              })();

              const barWidth = x.bandwidth();
              //const barHeight = innerHeight - y(obj['count']);
              // eslint-disable-next-line no-invalid-this
              // @ts-ignore
              // eslint-disable-next-line no-invalid-this
              const st = select(this);
              let that = st
                .append('rect');
              if (obj['selectedCount'] > 0) {
                st.append('line')
                  .style('stroke', 'orange')
                  .attr('x1', xPos - barWidth / 10)
                  .attr('y1', yPos)
                  .attr('x2', xPos - barWidth / 10)
                  .attr('y2', innerHeight - y(obj['selectedCount']) + yPos)
                  .style('stroke-width', barWidth / 7);
              }
              that = that.attr('class', 'bar')
                .attr('data-index', j)

                .attr('x', function(e: any) {
                  return x(d['name']);
                })
                .attr('width', x.bandwidth());
              // @ts-ignore
              (scope.dataFrame.currentRow[`${scope.aminoColumnNames[d['name'] - 1]}`] === obj['name']) ?
                that.style('stroke', 'black')
                  .style('stroke-width', Math.min(y(obj['count']), x.bandwidth()) / 10) : that.lower();

              that.style('fill', function(e: any) {
                return getColor(obj['name']);
              });

              if (innerHeight - y(obj['count']) > Math.min(x.bandwidth(), 10)) {
                st.append('text')
                  .text(function(e: any) {
                    return obj['name'];
                  })
                  .attr('x', function(e: any) {
                    return x(d['name']) + x.bandwidth() / 2;
                  })
                  .attr('y', function(e: any) {
                    let sum = 0;
                    arr.map((obj1: any, k: any) => {
                      if (k < j) {
                        sum = sum + obj1['count'];
                      }
                    });
                    return y(obj['count'] + sum) +
                                        Math.min(x.bandwidth(), 15) / 2 - 2 +
                                        (innerHeight - y(obj['count'])) / 2;
                  })
                  .attr('font-family', '"Open Sans", sans-serif')
                  .attr('font-size', `${Math.min(x.bandwidth(), 15)}px`)
                  .attr('fill', 'black')
                  .attr('text-anchor', 'middle');
              }
              that.attr('y', function(e: any) {
                let sum = 0;
                arr.map((obj1: any, k: any) => {
                  if (k < j) {
                    sum = sum + obj1['count'];
                  }
                });
                return y(obj['count'] + sum);
              })
                .attr('height', function(e: any) {
                  return innerHeight - y(obj['count']);
                })
                .on('mouseover', (event: any, e: any) => {
                  const newColor = getColor(obj['name']);
                  // @ts-ignore
                  that.style('fill', color(newColor).brighter([1]));
                  ui.tooltip.show(ui.divText(`${obj['name']}:${obj['count']}`), event.x, event.y);
                })
                .on('mouseout', () => {
                  const newColor = getColor(obj['name']);
                  that.style('fill', newColor);
                  ui.tooltip.hide();
                })
                .on('mousedown', (event: any, e: any) => {
                  // @ts-ignore
                  scope.dataFrame.selection.handleClick((i) => {
                    let aminoName = obj['name'];
                    //let selected = true;
                    if (obj['name'].endsWith('_unselected')) {
                      aminoName = obj['name'].substring(0, obj['name'].length - 11);
                      //selected = false;
                    }
                    // @ts-ignore
                    return (aminoName === (scope.dataFrame.getCol(`${scope.aminoColumnNames[d['name'] - 1]}`).get(i)));

                    //&& (scope.dataFrame.selection.get(i) === selected);
                  }, event);
                });
            });
          },
        );
    }

    onPropertyChanged(property: DG.Property) {
      super.onPropertyChanged(property);
      this.render();
    }

    register(args: DG.GridCellRenderArgs) {
      this.registered[args.cell.tableColumn!.name] = args.cell;
    }

    unregister(name: string) {
      if (this.registered[name]) {
        delete this.registered[name];
      }
    }


    renderBar(name: string) {
      if (!(this.registered[name]) || !(this.tableCanvas)) {
        return;
      }
      const cell = this.registered[name];
      const rect = cell.bounds;
      this.renderBarToCanvas(this.tableCanvas.getContext('2d')!, cell, rect.x, rect.y, rect.width, rect.height);
    }

    highlight(cell: DG.GridCell, offsetX:number, offsetY:number) {
      const colName = cell.tableColumn?.name;
      if (!colName) {
        return;
      }
      const margin = 0.2;
      const bound = cell.bounds;
      const x = bound.x + bound.width * margin;
      const y = 0 + 200 * margin / 4;
      const w = bound.width - bound.width * margin * 2;
      const h = 200 - 200 * margin / 2;
      const barData = this.barStats[colName];
      let sum = 0;
      barData.forEach((obj) => {
        sum += obj['count'];
      });
      let curSum = 0;

      barData.forEach((obj, index) => {
        const sBarHeight = h * obj['count'] / this.max;
        if (offsetX>=x &&
            offsetY>=y + h * (this.max - sum + curSum) / this.max &&
            offsetX<=x+w &&
            offsetY<=y + h * (this.max - sum + curSum) / this.max + sBarHeight) {
          this.highlighted = {'colName': colName, 'aaName': obj['name']};
          return;
        }

        curSum += obj['count'];
      });
      return;
    }

    unhighlight() {
      this.highlighted = null;
    }

    beginSelection(event:any) {
      if (!this.highlighted || !this.dataFrame) {
        return;
      }
      this.dataFrame!.selection.handleClick((i) => {
        //let selected = true;
        // @ts-ignore
        return this.highlighted!['aaName'] === (this.dataFrame.getCol(this.highlighted!['colName']).get(i));

        //&& (scope.dataFrame.selection.get(i) === selected);
      }, event);
    }
}
