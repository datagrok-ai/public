class StackedBarChart extends DG.JsViewer {
    constructor() {
        super();
        this.dataColumnPrefix = this.string('dataColumnPrefix', 'a');
        this.dataEmptyAA = this.string('dataEmptyAA', '-');
        this.valueAggrType = this.string('valueAggrType', 'avg', {choices: ['avg', 'count', 'sum']});
        this.initialized = false;
    }

    // Additional chart settings
    init() {
        let groups = [
            [['C', 'U'], 'yellow'],
            [['G', 'P'], 'red'],
            [['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'], 'all_green'],
            [['R', 'H', 'K'], 'light_blue'],
            [['D', 'E'], 'dark_blue'],
            [['S', 'T', 'N', 'Q'], 'orange']]
        this.ord = {}

        let i = 0
        for (const g in groups) {
            i++
            for (const obj in groups[g][0]) {
                this.ord[groups[g][0][obj]] = i
            }
        }

        this.margin = {top: 10, right: 10, bottom: 50, left: 10};
        this.yScale = scaleLinear();
        this.xScale = scaleBand();
        this.data = {};
        this.colors = {}
        this.palete = DG.Color.categoricalPalette.slice()
        this.palete.splice(this.palete.indexOf(4286545791), 1);
        this.colorScale = scaleOrdinal(this.palete);
        this.selection_mode = false
        this.aminoColumnNames = []
        let cp = ChemPalette.get_datagrock()
        this.getColor = (c = null) => {
            //return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
            return c in cp ? cp[c] : 'rgb(0,0,0)'
        }
        console.error(this.palete)
        this.initialized = true;
        for (let i = 0; i < 100; i++) {
            console.error(this.getColor(i))
        }
    }

    // Stream subscriptions
    onTableAttached() {
        this.init();
        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => this.render(false)));
        this.render();
    }

    // Cancel subscriptions when the viewer is detached
    detach() {
        this.subs.forEach(sub => sub.unsubscribe());
    }

    render(computeData = true) {
        let selection_mode = false
        if (computeData) {
            this.colors[this.dataEmptyAA] = "#f0f0f0"
            // Empty the data array
            this.data = {};
            this.aminoColumnNames = []
            this.dataFrame.columns.names().forEach((name) => {
                if (name.startsWith(this.dataColumnPrefix)) {
                    if (/^\d+$/.test(name.slice(this.dataColumnPrefix.length))) {
                        this.aminoColumnNames.push(name)
                    }
                }
            })

            this.aggregatedTables = {}
            this.aggregatedTablesUnselected = {}
            let buf1 = this.dataFrame.selection.getBuffer()
            let buf2 = this.dataFrame.filter.getBuffer()
            let resbuf = new Int32Array(this.dataFrame.rowCount)

            for (var i=0; i<buf2.length; i++) {
                resbuf[i] = buf1[i] & buf2[i];
            }



            let mask = DG.BitSet.fromBytes(resbuf.buffer,this.dataFrame.rowCount)
            if (mask.trueCount !== this.dataFrame.filter.trueCount) {
                this.selection_mode = true
                this.aminoColumnNames.forEach(name => {
                    this.aggregatedTables[name] = this.dataFrame
                        .groupBy([name])
                        .whereRowMask(this.dataFrame.filter)
                        .add('count', name, `${name}_count`)
                        .aggregate();
                    let buf1 =  this.dataFrame.selection.clone().invert().getBuffer()
                    let buf2 = this.dataFrame.filter.getBuffer()
                    let resbuf = new Int32Array(this.dataFrame.rowCount)

                    for (var i=0; i<buf2.length; i++) {
                        resbuf[i] = buf1[i] & buf2[i];
                    }



                    let mask = DG.BitSet.fromBytes(resbuf.buffer,this.dataFrame.rowCount)
                    this.aggregatedTablesUnselected[name] =  this.dataFrame
                        .groupBy([name])
                        .whereRowMask(mask)
                        .add('count', name, `${name}_count`)
                        .aggregate();
                });

            } else {
                this.selection_mode = false
                this.aminoColumnNames.forEach(name => {
                        this.aggregatedTables[name] = this.dataFrame
                            .groupBy([name])
                            .whereRowMask(this.dataFrame.filter)
                            .add('count', name, `${name}_count`)
                            .aggregate();
                    }
                )
            }

            this.data[`amino_count`] = {}
            this.data["colls"] = []
            let colobj = {}
            let empty_aa = null;


            for (const [name, df] of Object.entries(this.aggregatedTables)) {
                colobj = {"name": name.substring(1), "data": []}
                this.data["colls"].push(colobj)
                let unselected_row_index = 0

                for (let i = 0; i < df.rowCount; i++) {
                    let amino = df.getCol(name).get(i);
                    let amino_count = df.getCol(`${name}_count`).get(i)
                    if ((!amino) || amino in [this.dataEmptyAA]) {
                        //empty_aa = {"name": amino, "count": amino_count}
                        continue;
                    }
                    colobj["data"].push({"name": amino, "count": amino_count, "isSelected": this.selection_mode})
                    let unselected_amino_count = 0
                    let unselected_amino = ''

                    if ((name in this.aggregatedTablesUnselected) && (this.aggregatedTablesUnselected[name].rowCount > unselected_row_index)) {
                        unselected_amino = this.aggregatedTablesUnselected[name]
                            .getCol(name)
                            .get(unselected_row_index)

                        if (((!unselected_amino) || unselected_amino in [this.dataEmptyAA, ''])) {
                            unselected_row_index += 1
                            if (this.aggregatedTablesUnselected[name].rowCount > unselected_row_index) {
                                unselected_amino = this.aggregatedTablesUnselected[name]
                                    .getCol(name)
                                    .get(unselected_row_index)
                            }


                        }

                        if (unselected_amino === amino) {
                            unselected_amino_count = this.aggregatedTablesUnselected[name]
                                .getCol(`${name}_count`)
                                .get(unselected_row_index)
                            colobj["data"].slice(-1)[0]['count'] -= unselected_amino_count
                            colobj["data"].push({
                                "name": `${amino}`,
                                "count": unselected_amino_count,
                                "isSelected": false
                            })
                            unselected_row_index += 1
                        }
                    }

                    if (!(amino in this.data[`amino_count`])) {
                        this.colors[amino] = this.getColor(amino)
                        this.data[`amino_count`][amino] = 0
                    }

                    this.data[`amino_count`][amino]
                        += amino_count

                    if (amino === unselected_amino) {
                        this.data[`amino_count`][amino]
                            -= unselected_amino_count
                    }
                }
                colobj['data'] = colobj['data'].sort((o1, o2) => {
                    if (this.ord[o1['name']] > this.ord[o2['name']]) {
                        return -1
                    }
                    if (this.ord[o1['name']] < this.ord[o2['name']]) {
                        return 1
                    }

                    return 0
                })


            }

            if (empty_aa) {
                // colobj["data"].push(empty_aa)
                //empty_aa = null
            }

            this.data[`max`] = this.dataFrame.filter.trueCount
        }

        let width = this.root.parentElement.clientWidth;
        let height = this.root.parentElement.clientHeight;
        let innerWidth = width - this.margin.left - this.margin.right;
        let innerHeight = height - this.margin.top - this.margin.bottom;
        let colors = this.colors
        let scope = this

        $(this.root).empty();
        let svg = select(this.root).append("svg")
            .attr("width", width)
            .attr("height", height);
        let g = svg.append("g").attr("transform", `translate(${this.margin.left}, ${this.margin.top})`);

        const x = this.xScale
            .domain(this.data["colls"].map((d) => d['name']))
            .padding(0.3)
            .align(0.3)
            .rangeRound([0, innerWidth]);

        const y = this.yScale
            .domain([0, this.data["max"]]).nice()
            .rangeRound([innerHeight, 0]);

        // const aminoScale = scaleBand().domain(Object.entries(this.data[`amino_count`])
        //     .map((d) => (`${d[0]}:${d[1]}`)))
        //     .rangeRound([0, innerHeight]);

        let cAxis = axisBottom(x);
        cAxis.tickSize(0);
        let colAxis = g.append("g").call(cAxis)
            .attr("transform", `translate(0, ${innerHeight + 10})`);
        colAxis.attr("font-family", "common sans")
            .attr("fill", "black")
        colAxis.select(".domain").remove();


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

        g.selectAll('.group')
            .data(this.data["colls"])
            .attr('class', 'group')
            .enter().append('g')
            .each(function (d, i) {
                d['data'].map((obj, j, arr) => {
                    let that = select(this)
                        .append('rect')

                    that = that.attr('class', 'bar')
                        .attr('data-index', j)

                        .attr('x', function (e) {
                            return x(d['name']);
                        })
                        .attr('width', x.bandwidth());
                    that = obj['isSelected'] ?
                        that.style("stroke", "black")
                            .style("stroke-width", Math.min(y(obj['count']), x.bandwidth()) / 13) : that.lower();

                    that.style('fill', function (e) {
                        return (obj['name'] in colors) ? colors[obj['name']] : `rgb(127, 127, 127)`
                    })

                    if (innerHeight - y(obj['count']) > Math.min(x.bandwidth(), 10)) {
                        select(this).append("text")
                            .text(function (e) {
                                return obj['name'];
                            })
                            .attr("x", function (e) {
                                return x(d['name']) + x.bandwidth() / 2;
                            })
                            .attr("y", function (e) {
                                let sum = 0;
                                arr.map((obj_1, k) => {
                                    if (k < j) {
                                        sum = sum + obj_1['count'];
                                    }
                                });
                                return y(obj['count'] + sum) + Math.min(x.bandwidth(), 15) / 2 - 2 + (innerHeight - y(obj['count'])) / 2;
                            })
                            .attr("font-family", "common sans")
                            .attr("font-size", `${Math.min(x.bandwidth(), 15)}px`)
                            .attr("fill", "black")
                            .attr("text-anchor", "middle");
                    }
                    that.attr('y', function (e) {
                        let sum = 0;
                        arr.map((obj_1, k) => {
                            if (k < j) {
                                sum = sum + obj_1['count'];
                            }
                        });
                        return y(obj['count'] + sum);
                    })
                        .attr('height', function (e) {
                            return innerHeight - y(obj['count']);
                        })
                        .on('mouseover', (event, e) => {
                            let new_color = (obj['name'] in colors) ? colors[obj['name']] : `rgb(127, 127, 127)`
                            that.style('fill', color(new_color).brighter([1]));
                            ui.tooltip.show(ui.divText(`${obj['name']}:${obj['count']}`), event.x, event.y);

                        })
                        .on('mouseout', () => {
                            let new_color = (obj['name'] in colors) ? colors[obj['name']] : `rgb(127, 127, 127)`
                            that.style('fill', new_color);
                            ui.tooltip.hide()

                        })
                        .on('mousedown', (event, e) => {
                            scope.dataFrame.selection.handleClick(i => {
                                let amino_name = obj['name']
                                let selected = true;
                                if (obj['name'].endsWith('_unselected')) {
                                    amino_name = obj['name'].substring(0, obj['name'].length - 11)
                                    selected = false
                                }
                                let res = (amino_name === (scope.dataFrame.getCol(`a${d['name']}`).get(i)))
                                return res
                                //&& (scope.dataFrame.selection.get(i) === selected);
                            }, event);
                        });
                });
            });
    }

    onPropertyChanged(property) {
        super.onPropertyChanged(property);
        this.render();

    }
}
