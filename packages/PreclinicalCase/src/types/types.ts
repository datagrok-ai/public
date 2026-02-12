import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

export type StudyConfig = {
    protocol?: string;
    name?: string;
    description?: string;
    totalSubjects?: number;
    startDate?: dayjs.Dayjs;
    endDate?: dayjs.Dayjs;
    other?: {[key: string]: string};
    fieldsDefinitions?: {[key: string]: string};
}

export type TableView = {
    view: DG.TableView | DG.View,
    helper: any,
}
