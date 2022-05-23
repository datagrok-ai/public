import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class PersonView extends DG.ViewBase {

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {

        grok.data.query(`CDM:Person`)
            .then(personDf => {
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