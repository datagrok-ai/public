import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { createPivotedDataframe, filterNulls, getUniqueValues, getVisitNamesAndDays } from '../data-preparation/utils';
import { AE_START_DAY, AE_TERM, CON_MED_START_DAY, INV_DRUG_NAME, LAB_RES_N, LAB_TEST, SUBJECT_ID, VISIT_DAY, VISIT_NAME, VISIT_START_DATE, VS_RES_N, VS_TEST } from '../columns-constants';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';
import { createPropertyPanel } from '../panels/panels-service';


export class VisitsViewHelper {

    tv: DG.DataFrame;
    sv: DG.DataFrame;
    pivotedSv: DG.DataFrame;
    sortedVisitNamesAndDays: any;
    sortedVisitNames: any;
    patientVisit = new PatientVisit();
    studyVisit = new StudyVisit();
    name = 'Visits';
    totalVisits = {};
    proceduresAtVisit = { 'lb': LAB_RES_N, 'ex': INV_DRUG_NAME, 'vs': VS_RES_N };
    eventsSinceLastVisit = { 'ae': AE_START_DAY, 'cm': CON_MED_START_DAY };
    subjSet = new Set();

    constructor() {
        this.tv = study.domains.tv.clone();
        this.sv = study.domains.sv.clone();
        filterNulls(this.sv, VISIT_DAY);
        this.pivotedSv = createPivotedDataframe(this.sv, [SUBJECT_ID], VISIT_NAME, VISIT_START_DATE, []);
        this.pivotedSv.columns.names().forEach(col => {
            if (this.pivotedSv.getCol(col).name !== SUBJECT_ID) {
                this.pivotedSv.getCol(col).tags.format = 'yyyy-MM-dd';
            }
            this.pivotedSv.getCol(col).name = col.replace(` avg(${VISIT_START_DATE})`, '');
        });
        this.sortedVisitNamesAndDays = getVisitNamesAndDays(this.tv, true);
        this.sortedVisitNames = this.sortedVisitNamesAndDays.map(it => it.name);
        let missingCols = this.sortedVisitNames.filter(it => !this.pivotedSv.columns.names().includes(it));
        missingCols.forEach(it => {
            this.pivotedSv.columns.addNewDateTime(it).init((i) => null);
        })
        this.pivotedSv = this.pivotedSv.clone(null, [SUBJECT_ID].concat(this.sortedVisitNames));
        this.pivotedSv.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
                if (this.pivotedSv.currentCol.name !== SUBJECT_ID) {
                    if (this.pivotedSv.currentRowIdx === -1) {
                        let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(this.pivotedSv.currentCol.name);
                        this.studyVisit.createStudyVisitPanel(currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
                    } else {
                        let subjId = this.pivotedSv.get(SUBJECT_ID, this.pivotedSv.currentRowIdx);
                        let currentPatientVisit = this.totalVisits[this.pivotedSv.currentCol.name][subjId]
                        currentPatientVisit.updateSubjectVisit();
                        createPropertyPanel(currentPatientVisit);
                    }
                }

            }, 100);

        });
        this.createTotalVisits();

        Object.keys(this.proceduresAtVisit).concat(Object.keys(this.eventsSinceLastVisit)).forEach(domain => {
            let df = DG.DataFrame.create();
            df.columns.addNewString(SUBJECT_ID);
            this.pivotedSv.columns.names().forEach(col => {
                if(col !== SUBJECT_ID){
                    df.columns.addNewInt(col);
                }
            });
            Array.from(this.subjSet).forEach(id => {
                df.rows.addNew();
                df.getCol(SUBJECT_ID).set(df.rowCount, id);
                Object.keys(this.totalVisits).forEach(key => {
                    df.getCol(key).set(df.rowCount, this.totalVisits[key]);
                })
            });                
            this[domain] = df;
        })
    }


    getCurrentAndPreviousVisits(colName: string) {
        let visit = this.sortedVisitNamesAndDays.filter(it => it.name === colName)[0];
        let visitIdx = this.sortedVisitNamesAndDays.findIndex(it => it.name === colName);
        let previousVisit = visitIdx > 0 ? this.sortedVisitNamesAndDays[visitIdx - 1] : null;
        return { current: visit, previous: previousVisit };
    }

    createTotalVisits() {
        let countDfs = this.datasetsWithNumberProceduresAtVisit();
        this.createInitialTotalVisits();
        this.updateProceduresAtVisitCount(countDfs);
        this.updateEventsSinceLastVisitCount();
    }

    private datasetsWithNumberProceduresAtVisit() {
        let countDfs = {};
        Object.keys(this.proceduresAtVisit).forEach(domain => {
            if (study.domains[domain]) {
                countDfs[domain] = study.domains[domain]
                    .groupBy([SUBJECT_ID, VISIT_NAME])
                    .count(this.proceduresAtVisit[domain])
                    .aggregate();
            }
        });
        return countDfs;
    }

    private createInitialTotalVisits() {
        this.pivotedSv.columns.names().forEach(colName => {
            if (colName !== SUBJECT_ID) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(colName);
                this.totalVisits[colName] = {};

                this.pivotedSv.getCol(SUBJECT_ID).categories.forEach(subjId => {
                    this.totalVisits[colName][subjId] = new PatientVisit();
                    this.totalVisits[colName][subjId].currentSubjId = subjId;
                    this.totalVisits[colName][subjId].currentVisitDay = currentVisit.day;
                    this.totalVisits[colName][subjId].currentVisitName = currentVisit.name;
                    this.totalVisits[colName][subjId].previsousVisitDay = previousVisit ? previousVisit.day : null;
                });
            }

        });
    }

    private updateProceduresAtVisitCount(countDfs: any) {
        Object.keys(countDfs).forEach(domain => {
            for (let i = 0; i < countDfs[domain].rowCount; i++) {
                let visitName = countDfs[domain].get(VISIT_NAME, i);
                if (!this.sortedVisitNames.includes(visitName)) {
                    continue;
                }
                let subjId = countDfs[domain].get(SUBJECT_ID, i);
                this.subjSet.add(subjId);
                this.totalVisits[visitName][subjId].eventsCount[domain] = countDfs[domain].get(this.proceduresAtVisit[domain], i);
            }
        });
    }

    private updateEventsSinceLastVisitCount() {
        Object.keys(this.eventsSinceLastVisit).forEach(domain => {
            if (study.domains[domain]) {
                for (let i = 0; i < study.domains[domain].rowCount; i++) {
                    let startDay = study.domains[domain].get(this.eventsSinceLastVisit[domain], i);
                    let subjId = study.domains[domain].get(SUBJECT_ID, i);
                    this.subjSet.add(subjId);
                    for(let z = 0; z < this.sortedVisitNamesAndDays.length - 1; z++) {
                        if(z > 0) {
                            if (startDay > this.sortedVisitNamesAndDays[z].day && startDay < this.sortedVisitNamesAndDays[z+1].day){
                                let visitName = this.sortedVisitNamesAndDays[z+1].name;
                                if (this.totalVisits[visitName][subjId].eventsCount[domain]) {
                                    this.totalVisits[visitName][subjId].eventsCount[domain] += 1;
                                } else {
                                    this.totalVisits[visitName][subjId].eventsCount[domain] = 1;
                                }
                                break;
                            }
                        }
                    }
                }               
            }
        })
    }

}