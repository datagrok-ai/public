import * as DG from "datagrok-api/dg";
import { DataFrame } from "datagrok-api/dg";
import { study } from '../clinical-study';
import { SUBJECT_ID, VISIT_DAY, VISIT_NAME } from "../columns-constants";
import { dataframeContentToRow } from "../data-preparation/utils";
import { createPropertyPanel } from "../panels/panels-service";

export class PatientVisit {

    domainsNamesDict = {
        'Drug exposure': 'ex',
        'Laboratory': 'lb',
        'Vital signs': 'vs',
        'AEs since last visit': 'ae',
        'CONMEDs since last visit': 'cm'
    }
    atVisit = ['ex', 'lb', 'vs'];
    sinceLastVisit = ['ae', 'cm'];
    currentSubjId = '';
    currentVisitDay = null;
    currentVisitName = '';
    previsousVisitDay = null;
    name = 'Patient Visit';
    eventsCount = {};

    constructor() {
    }

    updateSubjectVisit() {
        this.atVisit.forEach(it => {
            if (study.domains[it]) {
                this[it] = study.domains[it]
                    .groupBy(study.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${VISIT_NAME} = ${this.currentVisitName}`)
                    .aggregate();
            }
        });

        this.sinceLastVisit.forEach(it => {
            if (study.domains[it] && this.previsousVisitDay) {
                this[it] = study.domains[it]
                    .groupBy(study.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${it.toUpperCase()}STDY > ${this.previsousVisitDay} and ${it.toUpperCase()}STDY <= ${this.currentVisitDay}`)
                    .aggregate();
                    if(this.currentSubjId === '01-701-1023' && this.currentVisitName === 'AMBUL ECG PLACEMENT') {
                        console.log(dataframeContentToRow(this[it]));
                    }
            } else {
                this[it] = null;
            }
        });

    }



}