import * as DG from "datagrok-api/dg";
import { DataFrame } from "datagrok-api/dg";
import { study } from '../clinical-study';
import { SUBJECT_ID } from "../columns-constants";
import { createPropertyPanel } from "../panels/panels-service";

export class AEBrowserHelper{

    domains = ['ae', 'ex', 'cm'];
    dyDomains = ['vs', 'lb'];
    domainsWithoutVisitDays = ['mh'];
    domainsToExclude = ['dm', 'tv', 'sv'];
    additionalDomains = [];
    selectedAdditionalDomains = [];
    aeToSelect: DG.DataFrame;
    daysPriorAe = 5;
    currentSubjId = '';
    currentAeDay: number;
    name = 'AE Browser';

    constructor(dataFrame: DG.DataFrame){
        let presentDomains = [];
        this.aeToSelect = dataFrame;
        study.domains.all().forEach(it => {
            if(!this.domainsToExclude.includes(it.name)){
                this[ it.name ] = study.domains[ it.name ].clone();
                presentDomains.push(it.name);
                if (!this.domains.includes(it.name)){
                    this.additionalDomains.push(it.name);
                }
            }
        });
        this.domains = this.domains.filter(it => presentDomains.includes(it));
    }

    updateDomains() {
        this.domains.concat(this.selectedAdditionalDomains).forEach(domain => {
            const condition = this.domainsWithoutVisitDays.includes(domain) ?
                `${SUBJECT_ID} = ${this.currentSubjId}` :
                this.dyDomains.includes(domain) ?
                    `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}DY < ${this.currentAeDay} and ${domain.toUpperCase()}DY > ${this.currentAeDay - this.daysPriorAe}` :
                    `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}STDY < ${this.currentAeDay} and ${domain.toUpperCase()}ENDY > ${this.currentAeDay - this.daysPriorAe}`;
            this[domain] = study.domains[domain]
                .clone()
                .groupBy(study.domains[domain].columns.names())
                .where(`${condition}`)
                .aggregate();
        })
    }

    createAEBrowserPanel(){
        this.updateDomains();
        createPropertyPanel(this);
    }

}