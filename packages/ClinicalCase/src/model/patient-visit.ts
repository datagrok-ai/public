import * as DG from "datagrok-api/dg";
import { DataFrame } from "datagrok-api/dg";
import { ClinicalDomains } from '../clinical-study';
import { SUBJECT_ID, VISIT_DAY, VISIT_NAME } from "../columns-constants";
import { dataframeContentToRow } from "../data-preparation/utils";

export class PatientVisit {

    domainsNamesDict = {
        'Drug exposure': 'ex',
        'Laboratory': 'lb',
        'Vital signs': 'vs',
        'AEs since last visit': 'ae',
        'CONMEDs since last visit': 'cm'
    };
    domains: ClinicalDomains;
    atVisit = ['ex', 'lb', 'vs'];
    sinceLastVisit = ['ae', 'cm'];
    currentSubjId = '';
    currentVisitDay = null;
    currentVisitName = '';
    previousVisitDay = null;
    name = 'Patient Visit';
    eventsCount = {};

    constructor(domains: ClinicalDomains) {
        this.domains = domains;
    }

    updateSubjectVisit(subjId: string, currentVisiDay: number, currentVisitName: string, previousVisitDay: number){
        this.currentSubjId = subjId;
        this.currentVisitDay = currentVisiDay;
        this.currentVisitName = currentVisitName;
        this.previousVisitDay = previousVisitDay
    }

    updateSubjectVisitDomains() {
        this.atVisit.forEach(it => {
            if (this.domains[it] && this.domains[it].columns.names().includes(VISIT_NAME)) {
                this[it] = this.domains[it]
                    .groupBy(this.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${VISIT_NAME} = ${this.currentVisitName}`)
                    .aggregate();
            }
        });

        this.sinceLastVisit.forEach(it => {
            if (this.domains[it] && this.previousVisitDay && this.domains[it].columns.names().includes(`${it.toUpperCase()}STDY`)) {
                this[it] = this.domains[it]
                    .groupBy(this.domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${it.toUpperCase()}STDY > ${this.previousVisitDay} and ${it.toUpperCase()}STDY <= ${this.currentVisitDay}`)
                    .aggregate();
            } else {
                this[it] = null;
            }
        });
    }

}