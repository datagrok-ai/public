import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {axisBottom, scaleBand, scaleLinear, select, color} from 'd3';
import {ChemPalette} from '../utils/chem-palette';
import $ from 'cash-dom';
import {GridCell, Property, Widget} from 'datagrok-api/dg';


export function addViewerToHeader(grid: DG.Grid, viewer: Promise<Widget>) {
  const wrapper = ui.div([]);
  let barchartRoot: HTMLElement = ui.div([]);
  let left = 200;
  let right = 200;
  let bottom = 300;
  console.error(wrapper);

  wrapper.style.width = '300px';
  wrapper.style.height = '200px';
  wrapper.style.position = 'absolute';
  wrapper.style.left = '200px';
  wrapper.setAttribute('position', 'absolute;');
  grid.setOptions({'colHeaderHeight': 200});
  grid.root.appendChild(wrapper);

  viewer.then((viewer) => {
    const barchart = viewer as StackedBarChart;
    barchartRoot = barchart.root;
    barchartRoot.setAttribute('class', 'ui-div');
    wrapper.appendChild(barchart.root);
    barchart.aminoColumnNames.forEach((name) => {
            grid.columns.byName(name)!.width = 20;
    });
    grid.onCellRender.subscribe( (args)=> {
      if (args.cell.isColHeader) {
        args.preventDefault();
      }
    });
    grid.onCellPrepare((cell: GridCell) => {
      if (cell.isColHeader && barchart.aminoColumnNames.length > 0) {
        if (cell.gridColumn.name == barchart.aminoColumnNames.at(-1)) {
          right = cell.bounds.right;
          bottom = cell.bounds.bottom;
          wrapper.style.width = `${left - right}px`;
          wrapper.style.height = `${bottom}px`;
          wrapper.style.left = `${left}px`;
          barchartRoot.style.width = `${left - right}px`;
          barchartRoot.style.height = `${bottom}px`;
          barchartRoot.style.left = `${left}px`;
        }
        if (cell.gridColumn.name == barchart.aminoColumnNames.at(0)) {
          left = cell.bounds.left;
          bottom = cell.bounds.bottom;
          wrapper.style.width = `${left - right}px`;
          wrapper.style.height = `${bottom}px`;
          wrapper.style.left = `${left}px`;
          barchartRoot.style.width = `${left - right}px`;
          barchartRoot.style.height = `${bottom}px`;
          barchartRoot.style.left = `${left}px`;
        }
      }
    });
  });
}


export class StackedBarChart extends DG.JsViewer {
    public dataColumnPrefix: string;
    public dataEmptyAA: string;
    public initialized: boolean;
    public valueAggrType: string;
    private ord: { [Key: string]: number; } = {};
    private margin: { top: number; left: number; bottom: number; right: number } = {
      top: 10,
      right: 10,
      bottom: 50,
      left: 10,
    };
    private yScale: any;
    private xScale: any;
    private data: { 'name': number, 'data': { 'name': string, 'count': number, 'selectedCount': number }[] }[] = [];
    private selectionMode: boolean = false;
    public aminoColumnNames: string[] = [];
    // @ts-ignore
    private getColor: ((c?: string) => string);
    private aminoColumnIndices: { [Key: string]: number; } = {};
    private aggregatedTables: { [Key: string]: DG.DataFrame; } = {};
    private aggregatedTablesUnselected: { [Key: string]: DG.DataFrame; } = {};
    private max = 0;

    constructor() {
      super();

      this.dataColumnPrefix = this.string('dataColumnPrefix', 'a');
      this.dataEmptyAA = this.string('dataEmptyAA', '-');
      this.valueAggrType = this.string('valueAggrType', 'avg', {choices: ['avg', 'count', 'sum']});
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
      const cp = ChemPalette.get_datagrok();
      this.getColor = (c = '') => {
        //return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
        return c in cp ? cp[c] : 'rgb(0,0,0)';
      };
    }

    // Stream subscriptions
    onTableAttached() {
      this.init();
      if (this.dataFrame) {
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.onCurrentRowChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
      }
      this.render();
    }

    // Cancel subscriptions when the viewer is detached
    detach() {
      this.subs.forEach((sub) => sub.unsubscribe());
    }

    render(computeData = true) {
      let df = DG.DataFrame.create();
      if (this.dataFrame) {
        df = this.dataFrame;
      } else {
        throw new Error('undefined dataframe');
      }

      if (computeData) {
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
        for (const [name, df] of Object.entries(this.aggregatedTables)) {
          const colObj: {
                    'name': number, 'data':
                        { 'name': string, 'count': number, 'selectedCount': number }[]
                } =
                    {'name': this.aminoColumnIndices[name], 'data': []};
          this.data.push(colObj);
          let unselectedRowIndex = 0;
          for (let i = 0; i < df.rowCount; i++) {
            const amino = df.getCol(name).get(i);
            const aminoCount = df.getCol(`${name}_count`).get(i);
            if ((!amino) || amino === this.dataEmptyAA) {
              continue;
            }
            const aminoObj = {'name': amino, 'count': aminoCount, 'selectedCount': 0};
            colObj['data'].push(aminoObj);


            if (name in this.aggregatedTablesUnselected) {
              if (amino != this.aggregatedTablesUnselected[name].getCol(name).get(unselectedRowIndex)) {
                unselectedRowIndex++;
              }
              aminoObj['selectedCount'] = this.aggregatedTablesUnselected[name]
                .getCol(`${name}_count`)
                .get(unselectedRowIndex);

              unselectedRowIndex++;
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
      }
      this.max = df.filter.trueCount;

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

      const getColor = this.getColor;
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

    onPropertyChanged(property: Property) {
      super.onPropertyChanged(property);
      this.render();
    }
}
