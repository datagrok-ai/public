import * as DG from "datagrok-api/dg";

export const VIEWS = [];

export const PERSON_ID = 'person_id';

export const conceptColumnNames = {
    conceptId: "Concept Id",
    conceptPath: "Name",
    numPersons: "Person Count",
    percentPersons: "Prevalence",
    recordsPerPerson: "Records Per Person"
}
export const PREVALENCE_TABLE_NAME = 'Prevalence';

export const CHARTS = [
    {
        name: 'PrevalenceByGenderAgeYear',
        chartName: 'Prevalence by Gender, Age, Year',
        type: DG.VIEWER.LINE_CHART,
        series: 'V',
        options: {
            split: 'trellis_name',
            xColumnName: 'x_calendar_year',
            yColumnNames: [ 'y_prevalence_1000_pp' ],
        }

    },
    {
        name: 'PrevalenceByMonth',
        chartName: 'Prevalence by Month',
        type: DG.VIEWER.LINE_CHART,
        options: {
            xColumnName: 'x_calendar_month',
            yColumnNames: [ 'y_prevalence_1000_pp' ],
        }
    },
    {
        name: 'ByType',
        chartName: 'Type',
        type: DG.VIEWER.PIE_CHART,
        options: {
            category: 'concept_name',
            segmentAngle: 'count_value',
            segmentAngleAggrType: 'avg',
        }
    }

]