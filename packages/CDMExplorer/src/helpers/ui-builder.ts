import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { CHARTS, conceptColumnNames } from '../constants';
import { _package } from '../package';
import { updateDivInnerHTML } from '../utils';
import { reportTitle, seriesTitle, viewerTitle, viewerTitleWithHeight } from '../styles';

export class UIBuilder {

    name: string;
    root: HTMLElement;
    charts: string[];

    constructor(name: string, viewRoot: HTMLElement, charts: string[]) {
        this.name = name;
        this.root = viewRoot;
        this.charts = charts;
        this.createUI();

    }

    private async createUI() {
        let viewContainer = ui.splitV([]);
        const nameWithoutSpaces = (this.name as any).replaceAll(' ', '');
        grok.data.query(`CDM:${nameWithoutSpaces}`).then(res => {
            Object.keys(conceptColumnNames).forEach(key => res.getCol(key).name = conceptColumnNames[ key ]);
            this.root.className = 'grok-view ui-box';
            viewContainer.append(res.plot.grid().root);
            this.root.appendChild(viewContainer);
            res.onCurrentRowChanged.subscribe(() => {
                Promise.all(this.charts.map(it =>
                    grok.data.query(`CDM:${it}For${nameWithoutSpaces}`,
                        { conceptId: res.get(conceptColumnNames[ 'conceptId' ], res.currentRowIdx) })))
                    .then(resultDfs => {
                        updateDivInnerHTML(viewContainer, res.plot.grid().root);
                        viewContainer.append(ui.divText(`${res.get(conceptColumnNames[ 'conceptPath' ], res.currentRowIdx)} Drilldown Report`, reportTitle))
                        resultDfs.forEach(df =>{
                            let pos = df.name.lastIndexOf(`For ${this.name}`);
                            let chartNname = df.name.substring(0,pos).replaceAll(' ', '');
                            let chartOptions = CHARTS.filter(chart => chart.name === chartNname)[0];
                            if(chartOptions.series){
                                let seriesContainer = chartOptions.series === 'H' ? ui.splitH([]) : null;
                                viewContainer.append(ui.divText(`${chartOptions.chartName}`, viewerTitleWithHeight))
                                let seriesuniqueValues = df.getCol('series_name').categories;
                                seriesuniqueValues.forEach(it => {
                                    let seriesDf = df.groupBy(df.columns.names())
                                    .where(`series_name = ${it}`)
                                    .aggregate();
                                    let seriesChart = DG.Viewer.fromType(chartOptions.type, seriesDf, chartOptions.options);
                                    if(seriesContainer){
                                        seriesContainer.append(seriesChart.root);
                                    } else {
                                        viewContainer.append(seriesChart.root);
                                    }
                                    seriesChart.root.prepend(ui.divText(`${it}`, seriesTitle));
                                })
                                if(seriesContainer){
                                    viewContainer.append(seriesContainer);
                                }
                            } else {
                                let chart = DG.Viewer.fromType(chartOptions.type, df, chartOptions.options);
                                chart.root.prepend(ui.divText(`${chartOptions.chartName}`, viewerTitle));
                                viewContainer.append(chart.root);
                            }
                        })
                    })
            });
        })
    }
}