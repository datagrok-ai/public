import * as DG from "datagrok-api/dg";
import { DataFrame } from "datagrok-api/dg";
import { study } from '../clinical-study';
import { SUBJECT_ID, VISIT_DAY, VISIT_NAME } from "../columns-constants";
import { createPropertyPanel } from "../panels/panels-service";

export class PatientVisit {

    domainsNamesDict = {
        'Drug exposure': 'ex',
        'Laboratory': 'lb',
        'Vitsl signs': 'vs',
        'New AEs since last visit': 'ae',
        'New CONMEDs since last visit': 'dm'
    }
    atVisit = ['ex', 'lb', 'vs'];
    sinceLastVisit = ['ae', 'cm'];
    currentSubjId = '';
    currentVisitDay = null;
    currentVisitName = '';
    previsousVisitDay = null;
    name = 'Patient Visit';

    constructor() {
    }

    updateSubjectVisit() {
        this.atVisit.forEach(it => {
            if (study.domains[it]) {
                this[it] = study.domains[it].clone()
                    .groupBy(study.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${VISIT_NAME} = ${this.currentVisitName}`)
                    .aggregate()
            }
        });

        this.sinceLastVisit.forEach(it => {
            if (study.domains[it] && this.previsousVisitDay) {
                this[it] = study.domains[it].clone()
                    .groupBy(study.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${it.toUpperCase()}STDY > ${this.previsousVisitDay} and ${it.toUpperCase()}STDY <= ${this.currentVisitDay}`)
                    .aggregate()
            } else {
                this[it] = null;
            }
        });

    }

    createPatientVisitPanel(subjId: string, visitDay: number, visitName: string, previsousVisitDay: number){
        this.currentSubjId = subjId;
        this.currentVisitDay = visitDay;
        this.currentVisitName = visitName;
        this.previsousVisitDay = previsousVisitDay;
        this.updateSubjectVisit();
        createPropertyPanel(this);
    }

}