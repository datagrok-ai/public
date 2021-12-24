import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { requiredColumnsByView } from '../constants';
import { createPivotedDataframe, filterNulls, getUniqueValues, getVisitNamesAndDays } from '../data-preparation/utils';
import { SUBJECT_ID, VISIT_DAY, VISIT_NAME, VISIT_START_DATE } from '../columns-constants';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';


export class VisitsView extends DG.ViewBase implements ILazyLoading {

    tv: DG.DataFrame;
    sv: DG.DataFrame;
    pivotedSv: DG.DataFrame;
    sortedVisitNamesAndDays: any;
    sortedVisitNames: any;
    patientVisit = new PatientVisit();
    studyVisit = new StudyVisit();

    constructor(name) {
        super({});
        this.name = name;
    }

    loaded: boolean;

    load(): void {
        checkMissingDomains(requiredColumnsByView[this.name], this);
    }

    createView(): void {
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
                if (this.pivotedSv.currentRowIdx === -1) {
                    let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits();
                    this.studyVisit.createStudyVisitPanel(currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
                } else {
                    let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits();
                    let subjId = this.pivotedSv.get(SUBJECT_ID, this.pivotedSv.currentRowIdx);
                    this.patientVisit.createPatientVisitPanel(subjId, currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
                }
            }, 100);

        });
        let grid = this.pivotedSv.plot.grid().root;

        this.root.className = 'grok-view ui-box';
        this.root.append(ui.splitV([
            grid
        ]));
    }

    private getCurrentAndPreviousVisits(){
        let visit = this.sortedVisitNamesAndDays.filter(it => it.name === this.pivotedSv.currentCol.name)[0];
        let visitIdx = this.sortedVisitNamesAndDays.findIndex(it => it.name === this.pivotedSv.currentCol.name);
        let previousVisit = visitIdx > 0 ? this.sortedVisitNamesAndDays[visitIdx - 1] : null;
        return {current: visit, previous: previousVisit};
    }


}