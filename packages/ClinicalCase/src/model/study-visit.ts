import * as DG from "datagrok-api/dg";
import { DataFrame, DateTime } from "datagrok-api/dg";
import { ClinicalDomains, ClinicalStudy} from '../clinical-study';
import { INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, INV_DRUG_NAME, SUBJECT_ID, VISIT_DAY, VISIT_NAME, VISIT_START_DATE } from "../columns-constants";
import { addColumnWithDrugPlusDosage } from "../data-preparation/data-preparation";

export class StudyVisit {

    domains: ClinicalDomains;
    currentVisitDay = null;
    currentVisitName = '';
    previsousVisitDay = null;
    name = 'Study Visit';
    totalPatients = 0;
    minVisitDate: string;
    maxVisitDate: string;
    visitDataframe: DG.DataFrame;
    aeSincePreviusVisit: DG.DataFrame;
    conmedSincePreviusVisit: DG.DataFrame;
    lbAtVisit: DG.DataFrame;
    exAtVisit: DG.DataFrame;
    vsAtVisit: DG.DataFrame;
    extrtWithDoseColName = 'EXTRT_WITH_DOSE';

    constructor(domains: ClinicalDomains) {
        this.domains = domains;
    }

    updateStudyVisit (visitDay: number, visitName: string, previsousVisitDay: number) {   
        this.currentVisitDay = visitDay;
        this.currentVisitName = visitName;
        this.previsousVisitDay = previsousVisitDay;    
        this.visitDataframe = this.domains.sv.groupBy(this.domains.sv.columns.names())
        .where(`${VISIT_NAME} = ${this.currentVisitName}`)
        .aggregate();
        this.totalPatients = this.visitDataframe.getCol(SUBJECT_ID).stats.uniqueCount;
        this.minVisitDate = new Date(this.visitDataframe.getCol(VISIT_START_DATE).stats.min * 1e-3).toLocaleDateString();
        this.maxVisitDate = new Date(this.visitDataframe.getCol(VISIT_START_DATE).stats.max * 1e-3).toLocaleDateString();
        this.aeSincePreviusVisit = this.createEventSincePeviousVisitDf('ae');
        this.conmedSincePreviusVisit = this.createEventSincePeviousVisitDf('cm');
        if (this.domains.ex && this.domains.ex.columns.names().includes(VISIT_NAME)) {
            let ex = this.domains.ex.clone();
            if (ex.columns.names().includes(INV_DRUG_NAME)) {
                if ([INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS].every(it => ex.columns.names().includes(it))) {
                    addColumnWithDrugPlusDosage(ex, INV_DRUG_NAME, INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, this.extrtWithDoseColName);
                } else {
                    ex.col(INV_DRUG_NAME).name = this.extrtWithDoseColName;
                }
            }
            this.exAtVisit = ex.groupBy(ex.columns.names())
            .where(`${VISIT_NAME} = ${this.currentVisitName}`)
            .aggregate();
        };
        if (this.domains.lb && this.domains.lb.columns.names().includes(VISIT_NAME)) {
            this.lbAtVisit = this.domains.lb.groupBy(this.domains.lb.columns.names())
            .where(`${VISIT_NAME} = ${this.currentVisitName}`)
            .aggregate();
        };
        if (this.domains.vs && this.domains.vs.columns.names().includes(VISIT_NAME)) {
            this.vsAtVisit = this.domains.vs.groupBy(this.domains.vs.columns.names())
            .where(`${VISIT_NAME} = ${this.currentVisitName}`)
            .aggregate();
        };

    }

    private createEventSincePeviousVisitDf(domain: string) {
        if (this.domains[domain] && this.previsousVisitDay) {
            return this.domains[domain].groupBy(this.domains[domain].columns.names())
                .where(`${domain.toUpperCase()}STDY > ${this.previsousVisitDay} and ${domain.toUpperCase()}STDY <= ${this.currentVisitDay}`)
                .aggregate();
        }
        return null;
    }

}