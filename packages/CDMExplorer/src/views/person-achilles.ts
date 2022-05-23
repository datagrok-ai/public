import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { viewerTitle } from '../styles';

export class PersonAchillesView extends DG.ViewBase {

    analysisIds = {'Gender': 2, 'Race': 4, 'Ethnicity': 5};

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {

        let pieChartsDiv = ui.splitH([]);
        let historgam;
        Promise.all([
            grok.data.query(`CDM:genderEthnicityRace`), 
            grok.data.query(`CDM:year`)
            ]).then(dataFrames => {
            dataFrames.forEach(item => {
                if(item.name === 'Year') {
                    historgam = item.plot.line({
                        x: 'year',
                        yColumnNames: ['countValue'],
                        chartTypes: ['Stacked Bar Chart']});
                    historgam.root.prepend(ui.divText(item.name, viewerTitle));
                } else {
                    Object.keys(this.analysisIds).forEach(key => {
                        let df = item
                        .groupBy(item.columns.names())
                        .where(`analysisId = ${this.analysisIds[key]}`)
                        .aggregate();
                        let pieChart = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, df, {
                            category: 'conceptName',
                            segmentAngle: 'countValue',
                            segmentAngleAggrType: 'avg',
                        });
                        pieChart.root.prepend(ui.divText(key, viewerTitle));
                        pieChartsDiv.append(pieChart.root);
                    })
                }
            });
            this.root.className = 'grok-view ui-box';
            this.root.appendChild(ui.splitV([
                pieChartsDiv,
                historgam.root
            ]));
          });
    }
    
}