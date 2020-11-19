import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import * as d3 from 'd3';


export class ChordViewer extends DG.JsViewer {

    constructor() {
        super();
        this.initialized = false;
        this.numColumns = [];
        this.strColumns = [];
        this.inpCol;
        this.outCol;
        this.data = [];
        this.conf = {};
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
            this.inpCol = this.strColumns[0];
            this.outCol = this.strColumns[1];
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

        this.data = this.inpCol.categories.map(cat => { return {
            id: cat.toLowerCase(),
            label: cat,
            len: 0,
            color: "#80b1d3"
        }});

        for (let i = 0; i < this.inpCol.length; i++) {

            let inpId = this.inpCol.get(i).toLowerCase();
            let outId = this.outCol.get(i).toLowerCase();
            this.data.find(obj => obj.id === inpId).len += 100;

            // TODO: Calculate `start` and `end` based on numColumns
            this.chords.push({
                source: {
                    id: inpId,
                    start: 1,
                    end: 12
                },
                target: {
                    id: outId,
                    start: 1,
                    end: 12
                }
            });
        }
        console.log(this.data);
        console.log(this.chords);
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

        circos.layout(this.data);
        circos.chords('beta-track', this.chords);
        circos.render();
        document.getElementById('chart')
                .children[0].children[0]
                .setAttribute('viewBox', `${-size/2} ${-size/2} ${size*2} ${size*2}`);
    }
}
