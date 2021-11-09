import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from '../clinical-study';
import { createPropertyPanel } from '../panels/panels-service';
import { SUBJECT_ID } from '../columns-constants';
import { ILazyLoading } from '../lazy-loading/lazy-loading';

let filters = [ 'USUBJID', 'AESEV', 'AEBODSYS', 'AESTDY' ]

export class AeBrowserView extends DG.ViewBase implements ILazyLoading {

    domains = [ 'ae', 'cm', 'ex' ];
    additionalDomains = [];
    selectedAdditionalDomains = [];
    aeToSelect: DG.DataFrame;
    daysPriorAe = 5;
    currentSubjId = '';
    currentAeDay: number;

    constructor(name) {
        super({});
        this.name = name;
        this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/ae_browser.md';
    }

    loaded: boolean;
    
    load(): void {
               study.domains.all().forEach(it => {
            if(it.name !== 'dm'){
                this[ it.name ] = study.domains[ it.name ].clone();
                if (!this.domains.includes(it.name)){
                    this.additionalDomains.push(it.name);
                }
            }
        });
        this.aeToSelect = study.domains.ae.clone();

        this.root.className = 'grok-view ui-box';

        this.root.append(ui.splitH([
            this.getFilters(),
            this.aeToSelect.plot.grid().root
        ]))

        this.subscribeToCurrentRow();
    }

    private getFilters() {
        let chart = DG.Viewer.fromType('Filters', this.aeToSelect, {
            'columnNames': filters,
            'showContextMenu': false,
        }).root;
        chart.style.overflowY = 'scroll';
        return chart
    }

    private subscribeToCurrentRow() {
        this.aeToSelect.onCurrentRowChanged.subscribe(() => {
            this.currentSubjId = this.aeToSelect.get(SUBJECT_ID, this.aeToSelect.currentRowIdx);
            this.currentAeDay = this.aeToSelect.get('AESTDY', this.aeToSelect.currentRowIdx);
            this.updateDomains();
            createPropertyPanel(this);
        })
    }

    updateDomains() {
        this.domains.concat(this.selectedAdditionalDomains).forEach(domain => {
                const condition = domain === 'lb' ?
                `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}DY < ${this.currentAeDay} and ${domain.toUpperCase()}DY > ${this.currentAeDay - this.daysPriorAe}` :
                `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}STDY < ${this.currentAeDay} and ${domain.toUpperCase()}ENDY > ${this.currentAeDay - this.daysPriorAe}`;
                this[ domain ] = study.domains[ domain ]
                    .clone()
                    .groupBy(study.domains[ domain ].columns.names())
                    .where(`${condition}`)
                    .aggregate();
        })
    }
}