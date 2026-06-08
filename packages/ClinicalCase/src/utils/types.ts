import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

export const SDTM = 'SDTM';

export type ClinStudyConfig = {
    protocol?: string;
    name?: string;
    description?: string;
    totalSubjects?: number;
    startDate?: dayjs.Dayjs;
    endDate?: dayjs.Dayjs;
    other?: {[key: string]: string};
    standard?: string;
    standardVersion?: string;
    fieldsDefinitions?: {[key: string]: string};
}

export type ClinCaseTableView = {
    view: DG.TableView | DG.View,
    helper: any,
}
