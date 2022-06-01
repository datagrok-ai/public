import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { convertColToString, joinCohorts } from '../preprocessing/data-preparation';
import { cohorts } from '../cohorts';
import { PERSON_ID } from '../constants';
import { CDMViewBase } from '../model/cdmViewBase';

export class DemographicsView extends CDMViewBase {

    constructor(name) {
        super({});
        this.name = name;
    }

    async createView(): Promise<void> {
        ui.setUpdateIndicator(this.root, true);
        grok.data.query(`CDM:Person`)
            .then(personDf => {
                convertColToString(personDf, PERSON_ID);
                joinCohorts(personDf);
                grok.data.linkTables(cohorts.cohortsPivoted, personDf,
                    [ PERSON_ID ], [ PERSON_ID ],
                    [ DG.SYNC_TYPE.FILTER_TO_FILTER ]);
                let genderPieChart = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, personDf, {
                    category: 'gender',
                });
                let racePieChart = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, personDf, {
                    category: 'race',
                });
                let ethnicityPieChart = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, personDf, {
                    category: 'ethnicity',
                });
                let yearOfBirthHistogram = DG.Viewer.fromType(DG.VIEWER.HISTOGRAM, personDf, {
                    value: 'year_of_birth'
                });
                ui.setUpdateIndicator(this.root, false);
                this.root.append(ui.splitV([
                    ui.splitH([
                        genderPieChart.root,
                        racePieChart.root
                    ]),
                    ui.splitH([
                        ethnicityPieChart.root,
                        yearOfBirthHistogram.root
                    ]),
                    personDf.plot.grid().root
                ]));

            });
    }

}