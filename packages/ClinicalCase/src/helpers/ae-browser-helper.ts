import * as DG from "datagrok-api/dg";
import { DataFrame } from "datagrok-api/dg";
import { study } from '../clinical-study';
import { SUBJECT_ID } from "../constants";
import { createPropertyPanel } from "../panels/panels-service";

export class AEBrowserHelper{

    domains = [ 'ae', 'cm', 'ex' ];
    additionalDomains = [];
    selectedAdditionalDomains = [];
    aeToSelect: DG.DataFrame;
    daysPriorAe = 5;
    currentSubjId = '';
    currentAeDay: number;
    name = 'AE Browser';

    constructor(dataFrame: DG.DataFrame){
        this.aeToSelect = dataFrame;
        study.domains.all().forEach(it => {
            if(it.name !== 'dm'){
                this[ it.name ] = study.domains[ it.name ].clone();
                if (!this.domains.includes(it.name)){
                    this.additionalDomains.push(it.name);
                }
            }
        });
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

    createAEBrowserPanel(){
        this.updateDomains();
        createPropertyPanel(this);
    }

}