class StackedBarChart extends DG.JsViewer {
    constructor() {
        super();
        this.dataColumnPrefix = this.string('dataColumnPrefix', 'a');
        this.dataEmptyAA = this.string('dataEmptyAA', '-');
        this.initialized = false;
    }

    init() {
        this.margin = {top: 40, right: 10, bottom: 40, left: 10};
        this.yScale = scaleLinear();
        this.xScale = scaleBand();
        this.data = {};
        this.colors = {}
        this.palete = DG.Color.categoricalPalette.slice()
        this.palete.splice(this.palete.indexOf(4286545791), 1);
        this.colorScale = scaleOrdinal(this.palete);
        this.getColor = (c = null) => {
            return c ? DG.Color.toRgb(this.colorScale(c)) : 'rgb(127,127,127)'
        }
        this.initialized = true;
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
        if (computeData) {
            this.colors[this.dataEmptyAA] = "#f0f0f0"
            // Empty the data array
            this.data = {};
            let aminoColumnNames = []
            this.dataFrame.columns.names().forEach((name) => {
                if (name.startsWith(this.dataColumnPrefix)) {
                    if (/^\d+$/.test(name.slice(this.dataColumnPrefix.length))) {
                        aminoColumnNames.push(name)
                    }
                }
            })

            this.aggregatedTables = {}
            this.aggregatedTablesUnselected = {}
            if (this.dataFrame.selection.trueCount > 0) {
                aminoColumnNames.forEach(name => {
                    this.aggregatedTables[name] = this.dataFrame
                        .groupBy([name])
                        .whereRowMask(this.dataFrame.filter)
                        .add('count', name, `${name}_count`)
                        .aggregate();

                    this.aggregatedTablesUnselected[name] =
                        this.dataFrame.groupBy([name])
                            .whereRowMask(this.dataFrame.filter)
                            .whereRowMask(this.dataFrame.selection.invert())
                            .add('count', name, `${name}_count`)
                            .aggregate();
                    this.dataFrame.selection.invert()
                });

            } else {
                aminoColumnNames.forEach(name => {
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
                colobj = {"name": name, "data": []}
                this.data["colls"].push(colobj)
                let unselected_row_index = 0

                for (let i = 0; i < df.rowCount; i++) {
                    let amino = df.getCol(name).get(i);
                    let amino_count = df.getCol(`${name}_count`).get(i)
                    if (amino === this.dataEmptyAA) {
                        //empty_aa = {"name": amino, "count": amino_count}
                        continue;
                    }
                    colobj["data"].push({"name": amino, "count": amino_count})
                    let unselected_amino_count = 0
                    let unselected_amino = ''

                    if ((name in this.aggregatedTablesUnselected) && (df.rowCount > unselected_row_index)) {
                        unselected_amino = this.aggregatedTablesUnselected[name]
                            .getCol(name)
                            .get(unselected_row_index)

                        if (unselected_amino === this.dataEmptyAA && (df.rowCount > unselected_row_index + 1)) {
                            unselected_row_index += 1
                            unselected_amino = this.aggregatedTablesUnselected[name]
                                .getCol(name)
                                .get(unselected_row_index)
                        }

                        if (unselected_amino === amino) {
                            unselected_amino_count = this.aggregatedTablesUnselected[name]
                                .getCol(`${name}_count`)
                                .get(unselected_row_index)
                            colobj["data"].slice(-1)[0]['count'] -= unselected_amino_count
                            colobj["data"].push({"name": `${amino}_unselected`, "count": unselected_amino_count})
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

        const aminoScale = scaleBand().domain(Object.entries(this.data[`amino_count`])
            .map((d) => (`${d[0]}:${d[1]}`)))
            .rangeRound([0, innerWidth]);


        let colAxis = g.append("g").call(axisTop(x))
        colAxis.select(".domain").remove();

        let aminoAxis = g.append("g").call(axisBottom(aminoScale))
            .attr("transform", `translate(0, ${innerHeight})`);
        aminoAxis.select(".domain").remove();
        let colors = this.colors
        let scope = this

        g.selectAll('.group')
            .data(this.data["colls"])
            .attr('class', 'group')
            .enter().append('g')
            .each(function (d, i) {
                d['data'].map((obj, j, arr) => {
                    let that = select(this)
                        .append('rect')

                    that.attr('class', 'bar')
                        .attr('data-index', j)
                        .attr('x', function (e) {
                            return x(d['name']);
                        })
                        .attr('width', x.bandwidth())
                        .style('fill', function (e) {
                            return (obj['name'] in colors) ? colors[obj['name']] : `rgb(127, 127, 127)`
                        })
                        .attr('y', function (e) {
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
                                let res = (amino_name === (scope.dataFrame.getCol(d['name']).get(i)))
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
