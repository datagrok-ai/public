/* eslint-disable no-unused-vars */
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

export type ClinStudyConfig = {
    protocol?: string;
    name?: string;
    description?: string;
    totalSubjects?: number;
    startDate?: dayjs.Dayjs;
    endDate?: dayjs.Dayjs;
    other?: {[key: string]: string};
    standard?: CDISC_STANDARD;
    fieldsDefinitions?: {[key: string]: string};
}

export enum CDISC_STANDARD {
    SDTM = 'SDTM',
    SEND = 'SEND',
}

export type ClinCaseTableView = {
    view: DG.TableView | DG.View,
    helper: any,
}
