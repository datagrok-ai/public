import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { updateDivInnerHTML } from '../utils';
import { convertColToString, joinCohorts } from '../preprocessing/data-preparation';
import { PERSON_ID } from '../constants';
import $ from "cash-dom";

export class TimelinesView extends DG.ViewBase {

    links = {
        'Condition Occurrence': { key: PERSON_ID, start: 'condition_start_date', end: 'condition_end_date', event: 'condition' },
        'Drug Exposure': { key: PERSON_ID, start: 'drug_exposure_start_date', end: 'drug_exposure_end_date', event: 'drug' },
    }

    options = {
        splitByColumnName: 'key',
        startColumnName: 'start',
        endColumnName: 'end',
        colorByColumnName: 'domain',
        eventColumnName: 'event',
        showEventInTooltip: true
    }


    dataframes: DG.DataFrame[];
    selectedOptions: string[] = [];
    selectedDataframes: any;
    timelinesDiv = ui.box();
    filtersDiv = ui.box();
    resultTables: DG.DataFrame;
    filterColumns = [];
    filters: any;

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {
        Promise.all([
            grok.data.query(`CDM:conditionOccurrence`),
            grok.data.query(`CDM:drugExposure`)
        ]).then(dataframes => {
            this.dataframes = dataframes;
            this.dataframes.forEach(df => {
                convertColToString(df, PERSON_ID);
                joinCohorts(df);
            })
            this.selectedOptions = ['Drug Exposure'];
            let multiChoiceOptions = ui.multiChoiceInput('', this.selectedOptions as any, ['Drug Exposure', 'Condition Occurrence']);
            multiChoiceOptions.onChanged((v) => {
                this.selectedOptions = multiChoiceOptions.value;
                this.updateTimelinesPlot();
            });
            let customTitle = {
                style: {
                    'color': 'var(--grey-6)',
                    'margin-top': '8px',
                    'font-size': '16px',
                }
            };

            let viewerTitle = {
                style: {
                    'color': 'var(--grey-6)',
                    'margin': '12px 0px 6px 12px',
                    'font-size': '16px',
                }
            };

            this.root.className = 'grok-view ui-box';
            this.append(ui.splitH([
                ui.splitV([
                    ui.box(ui.panel([
                        ui.divText('Events', customTitle),
                        multiChoiceOptions.root
                    ]), { style: { maxHeight: '130px' } }),
                    ui.box(
                        ui.divText('Filters', viewerTitle), { style: { maxHeight: '40px' } }),
                    this.filtersDiv
                ], { style: { maxWidth: '250px', } }),
                this.timelinesDiv
            ]))

            this.updateTimelinesPlot();

        })
    }


    private prepare(domain: DG.DataFrame) {
        let info = this.links[domain.name];
        let df = domain;
        let t = df.clone(null, Object.keys(info).map(e => info[e]));
        let filterCols = df.columns.names().filter(it => !Object.values(info).includes(it));
        filterCols.forEach(key => { t.columns.addNewString(key).init((i) => df.get(key, i).toString()); })
        t.columns.addNew('domain', DG.TYPE.STRING).init(domain.name);
        t.columns.addNewFloat('rowNum').init((i) => i);
        for (let name in info)
            t.col(info[name]).name = name;
        return t;
    }

    private updateTimelinesTables() {
        this.resultTables = null;
        for (let dt of this.dataframes.filter(it => this.selectedOptions.includes(it.name))) {
            let t = this.prepare(dt);
            if (this.resultTables == null)
                this.resultTables = t;
            else
                this.resultTables.append(t, true);
        }
    }

    private updateTimelinesPlot() {
        this.updateTimelinesTables();
        if (this.resultTables) {
            this.filterColumns = this.resultTables.columns.names().filter(it => !['domain', 'rowNum', 'key', 'start', 'end', 'event'].includes(it))
            this.resultTables.plot.fromType(DG.VIEWER.TIMELINES, {
                paramOptions: JSON.stringify(this.options),
            }).then((v: any) => {
                v.setOptions({
                    splitByColumnName: 'key',
                    startColumnName: 'start',
                    endColumnName: 'end',
                    colorByColumnName: 'domain',
                    eventColumnName: 'event',
                    showEventInTooltip: true
                });
                $(v.root).css('position', 'relative')
                v.zoomState = [[0, 10], [0, 10], [90, 100], [90, 100]];
                v.render();
                this.updateTimelinesDivs(v.root, this.getFilters());
            });
        } else {
            this.updateTimelinesDivs('', '');
        }

    }


    private getFilters() {
        let chart = DG.Viewer.fromType('Filters', this.resultTables, {
            'columnNames': this.filterColumns,
            'showContextMenu': true,
        }).root;
        chart.style.overflowY = 'scroll';
        return chart
    }

    private updateTimelinesDivs(timelinesContent: any, filtersContent: any) {
        updateDivInnerHTML(this.timelinesDiv, timelinesContent);
        updateDivInnerHTML(this.filtersDiv, filtersContent);
    }


}