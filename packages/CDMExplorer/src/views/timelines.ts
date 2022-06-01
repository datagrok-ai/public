import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { updateDivInnerHTML } from '../utils';
import { convertColToString, joinCohorts } from '../preprocessing/data-preparation';
import { PERSON_ID } from '../constants';
import $ from "cash-dom";
import { cohorts } from '../cohorts';
import { CDMViewBase } from '../model/cdmViewBase';

export class TimelinesView extends CDMViewBase {

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
    resultTables: DG.DataFrame;
    domainsMultichoice: DG.InputBase;

    constructor(name) {
        super({});
        this.name = name;
    }

    async createView(): Promise<void> {
        ui.setUpdateIndicator(this.root, true);
        Promise.all([
            grok.data.query(`CDM:conditionOccurrence`),
            grok.data.query(`CDM:drugExposure`)
        ]).then(dataframes => {
            this.dataframes = dataframes;
            this.dataframes.forEach(df => {
                convertColToString(df, PERSON_ID);
            });
            this.selectedOptions = ['Drug Exposure'];
            this.domainsMultichoice = ui.multiChoiceInput('', this.selectedOptions as any, ['Drug Exposure', 'Condition Occurrence']);
            this.domainsMultichoice.onChanged((v) => {
                this.selectedOptions = this.domainsMultichoice.value;
                this.updateTimelinesPlot();
            });

            this.root.className = 'grok-view ui-box';
            this.append(ui.splitH([
                this.timelinesDiv
            ]))
            ui.setUpdateIndicator(this.root, false);
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
        ui.setUpdateIndicator(this.timelinesDiv, true);
        this.updateTimelinesTables();
        if (this.resultTables) {
            grok.data.linkTables(cohorts.cohortsPivoted, this.resultTables,
                [ PERSON_ID ], [ 'key' ],
                [ DG.SYNC_TYPE.FILTER_TO_FILTER ]);
            //grok.shell.addTableView(this.resultTables);
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
                this.updateTimelinesDivs(v.root);
            });
        } else {
            this.updateTimelinesDivs('');
        }

    }


    private updateTimelinesDivs(timelinesContent: any) {
        ui.setUpdateIndicator(this.timelinesDiv, false);
        updateDivInnerHTML(this.timelinesDiv, timelinesContent);
    }

    override async propertyPanel() {

        let acc = this.createAccWithTitle(this.name);
        acc.addPane('Domain', () => this.domainsMultichoice.root, true);

        return acc.root;
    }


}