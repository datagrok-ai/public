/* eslint-disable no-unused-vars */
import * as DG from 'datagrok-api/dg';

export type ClinStudyConfig = {
    protocol?: string;
    name?: string;
    description?: string;
    totalSubjects?: number;
    startDate?: string;
    endDate?: string;
    other?: {[key: string]: string};
    standard?: CDISC_STANDARD;
}

export enum CDISC_STANDARD {
    SDTM = 'SDTM',
    SEND = 'SEND',
}

export type ClinCaseTableView = {
    view: DG.TableView | DG.View,
    helper: any,
}
