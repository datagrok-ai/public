import { SUBJECT_ID, VISIT_NAME } from "../constants/columns-constants";

export class PatientVisit {

    domainsNamesDict = {
        'Drug exposure': 'ex',
        'Laboratory': 'lb',
        'Vital signs': 'vs',
        'AEs since last visit': 'ae',
        'CONMEDs since last visit': 'cm'
    };
    atVisit = ['ex', 'lb', 'vs'];
    sinceLastVisit = ['ae', 'cm'];
    currentSubjId = '';
    currentVisitDay = null;
    currentVisitName = '';
    previousVisitDay = null;
    name = 'Patient Visit';
    eventsCount = {};

    constructor() {
    }

    updateSubjectVisit(subjId: string, currentVisiDay: number, currentVisitName: string, previousVisitDay: number){
        this.currentSubjId = subjId;
        this.currentVisitDay = currentVisiDay;
        this.currentVisitName = currentVisitName;
        this.previousVisitDay = previousVisitDay
    }

    updateSubjectVisitDomains(domains: any) {
        this.atVisit.forEach(it => {
            if (domains[it] && domains[it].columns.names().includes(VISIT_NAME)) {
                this[it] = domains[it]
                    .groupBy(domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${VISIT_NAME} = ${this.currentVisitName}`)
                    .aggregate();
            }
        });

        this.sinceLastVisit.forEach(it => {
            if (domains[it] && this.previousVisitDay && domains[it].columns.names().includes(`${it.toUpperCase()}STDY`)) {
                this[it] = domains[it]
                    .groupBy(domains[it].columns.names())
                    .where(`${SUBJECT_ID} = ${this.currentSubjId} and ${it.toUpperCase()}STDY > ${this.previousVisitDay} and ${it.toUpperCase()}STDY <= ${this.currentVisitDay}`)
                    .aggregate();
            } else {
                this[it] = null;
            }
        });
    }

}