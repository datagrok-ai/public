import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import * as d3 from 'd3';
import { layoutConf } from './configuration.js';


export class ChordViewer extends DG.JsViewer {

    constructor() {
        super();

        // properties
        this.fromColumnName = this.string('fromColumnName');
        this.toColumnName = this.string('toColumnName');

        this.initialized = false;
        this.numColumns = [];
        this.strColumns = [];
        this.inpCol;
        this.outCol;
        this.data = [];
        this.conf = layoutConf;
        this.chords = [];
        this.chordConf = {};
    }

    init() {
        this.initialized = true;
    }

    testColumns() {
        return (this.strColumns.length >= 2 && this.numColumns.length >= 1);
    }

    onTableAttached() {
        this.init();

        let columns = this.dataFrame.columns.toList();
        this.strColumns = columns.filter(col => col.type === 'string');
        this.numColumns = columns.filter(col => ['double', 'int'].includes(col.type));

        if (this.testColumns()) {
            // TODO: Choose the most relevant columns
            this.fromColumnName = this.strColumns[0].name;
            this.toColumnName = this.strColumns[1].name;
            this.generateData();
        }

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.onSizeChanged, 50).subscribe((_) => this.render()));

        this.render();
    }

    onPropertyChanged(property) {
        super.onPropertyChanged(property);
        if (this.initialized && this.testColumns()) {
            this.render();
        }
    }

    detach() {
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    generateData() {

        this.aggregatedTable = this.dataFrame
            .groupBy([this.fromColumnName, this.toColumnName])
            .count('count')
            .aggregate();

        let fromCol = this.aggregatedTable.columns.byName(this.fromColumnName);
        let toCol = this.aggregatedTable.columns.byName(this.toColumnName);

        function toCircos(s) { return {
            id: s,
            label: s,
            len: 0,
            color: "#80b1d3"
        }}

        this.data = fromCol.categories.map(toCircos).concat(toCol.categories.map(toCircos));

        for (let i = 0; i < this.aggregatedTable.rowCount; i++) {
            // TODO: Calculate `start` and `end` based on numColumns
            this.chords.push({
                source: {
                    id: fromCol.get(i),
                    start: 1,
                    end: 12
                },
                target: {
                    id: toCol.get(i),
                    start: 1,
                    end: 12
                }
            });
        }
    }

    render() {

        if (!this.testColumns()) {
            this.root.innerText = "Not enough data to produce the result.";
            return;
        }

        this.root.classList.add('viewer-window');
        let size = Math.min(this.root.clientWidth, this.root.clientHeight);
        let svg = d3.select(this.root)
                    .append("div")
                        .attr('id', 'chart')
                        .attr('class', 'chord-diagram-container')
                        .attr('width', size)
                        .attr('height', size);

        let circos = Circos({
            container: '#chart',
            width: size,
            height: size
        });

        circos.layout(this.data, this.conf);
        circos.chords('beta-track', this.chords);
        circos.render();
        document.getElementById('chart')
                .children[0].children[0]
                .setAttribute('viewBox', `${-size/2} ${-size/2} ${size*2} ${size*2}`);
    }
}
