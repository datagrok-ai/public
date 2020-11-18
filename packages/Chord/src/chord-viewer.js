import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Circos from 'circos';
import * as d3 from 'd3';


export class ChordViewer extends DG.JsViewer {

    constructor() {
        super();
        this.initialized = false;
    }

    init() {
        this.initialized = true;
    }

    onTableAttached() {
        this.init();

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.onSizeChanged, 50).subscribe((_) => this.render()));

        this.render();
    }

    onPropertyChanged(property) {
        if (this.initialized)
            this.render();
    }

    detach() {
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    render() {

        // Make a new method #1 [Q: do column methods like isNumeric/isAlpha exist?]
        let columns = this.dataFrame.columns.toList();
        let strColumns = columns.filter(col => col.type === 'string');
        let numColumns = columns.filter(col => ['double', 'int'].includes(col.type));
        if (strColumns.length < 2 || numColumns.length < 1) {
            // TODO: Disable the viewer's icon in the menu
            this.root.innerText = "Not enough data to produce the result."
            return;
        }

        // Make a new method #2
        // TODO: Choose the appropriate columns
        let inpCol = strColumns[0];
        let outCol = strColumns[1];
        let numCol = numColumns[0];

        let data = [];
        let chords = [];

        for (let i = 0; i < inpCol.length; i++) {
            let label = inpCol.get(i);
            let inpId = label.toLowerCase();
            let outId = outCol.get(i).toLowerCase();
            // TODO: Calculate the length of a segment as inpCol.categories().length
            // (and update it with every next segment). 

            data.push({
                id: inpId,
                label: label,
                len: numCol.get(i),
                color: "#80b1d3"
            });

            // TODO: Calculate `start` and `end`
            chords.push({
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
        console.log(data);
        console.log(chords);


        this.root.style += "width: 100%; height: 100%;";
        let size = Math.min(this.root.clientWidth, this.root.clientHeight);
        let svg = d3.select(this.root)
                    .append("div")
                        .attr('id', 'chart')
                        .attr('class', 'chord-diagram-container')
                        .attr('width', size)
                        .attr('height', size)

        let circos = Circos({
            container: '#chart',
            width: size,
            height: size
        })

        circos.layout(data);
        circos.chords('beta-track', chords);
        circos.render();
        document.getElementById('chart')
                .children[0].children[0]
                .setAttribute('viewBox', `${-size/2} ${-size/2} ${size*2} ${size*2}`)
    }
}
