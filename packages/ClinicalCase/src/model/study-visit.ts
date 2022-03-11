import * as DG from "datagrok-api/dg";
import { INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, INV_DRUG_NAME, SUBJECT_ID, VISIT_NAME, VISIT_START_DATE } from "../constants/columns-constants";
import { addColumnWithDrugPlusDosage } from "../data-preparation/data-preparation";

export class StudyVisit {

    day = null;
    name = '';
    num: number;
    previsousVisitDay = null;
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

    constructor(name?: string, day?: string) {
        this.name = name;
        this.day = day;
    }

    updateStudyVisit (domains: any, visitDay: number, visitName: string, previsousVisitDay: number) {   
        this.day = visitDay;
        this.name = visitName;
        this.previsousVisitDay = previsousVisitDay;    
        this.visitDataframe = domains.sv.groupBy(domains.sv.columns.names())
        .where(`${VISIT_NAME} = ${this.name}`)
        .aggregate();
        this.totalPatients = this.visitDataframe.getCol(SUBJECT_ID).stats.uniqueCount;
        this.minVisitDate = new Date(this.visitDataframe.getCol(VISIT_START_DATE).stats.min * 1e-3).toLocaleDateString();
        this.maxVisitDate = new Date(this.visitDataframe.getCol(VISIT_START_DATE).stats.max * 1e-3).toLocaleDateString();
        this.aeSincePreviusVisit = this.createEventSincePeviousVisitDf(domains['ae']);
        this.conmedSincePreviusVisit = this.createEventSincePeviousVisitDf(domains['cm']);
        if (domains.ex && domains.ex.columns.names().includes(VISIT_NAME)) {
            let ex = domains.ex.clone();
            if (ex.columns.names().includes(INV_DRUG_NAME)) {
                addColumnWithDrugPlusDosage(ex, INV_DRUG_NAME, INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, this.extrtWithDoseColName);
            }
            this.exAtVisit = ex.groupBy(ex.columns.names())
            .where(`${VISIT_NAME} = ${this.name}`)
            .aggregate();
        };
        if (domains.lb && domains.lb.columns.names().includes(VISIT_NAME)) {
            this.lbAtVisit = domains.lb.groupBy(domains.lb.columns.names())
            .where(`${VISIT_NAME} = ${this.name}`)
            .aggregate();
        };
        if (domains.vs && domains.vs.columns.names().includes(VISIT_NAME)) {
            this.vsAtVisit = domains.vs.groupBy(domains.vs.columns.names())
            .where(`${VISIT_NAME} = ${this.name}`)
            .aggregate();
        };

    }

    private createEventSincePeviousVisitDf(domain: any) {
        const startDayColName = `${domain.name.toUpperCase()}STDY`;
        if (domain && this.previsousVisitDay && domain.col(startDayColName)) {
            return domain.groupBy(domain.columns.names())
                .where(`${startDayColName} > ${this.previsousVisitDay} and ${startDayColName} <= ${this.day}`)
                .aggregate();
        }
        return null;
    }

}