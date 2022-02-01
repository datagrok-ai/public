import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Viewers', () => {
    let v: DG.TableView;

    before(async () => {
       v = grok.shell.addTableView(grok.data.demo.demog(5000));
       v.barChart();
       v.scatterPlot();
       v.histogram();
       v.pieChart();
       v.lineChart();
       v.trellisPlot();
       v.matrixPlot();
       v.scatterPlot3d();
       v.densityPlot();
       v.pcPlot();
       v.wordCloud();
       v.networkDiagram();
       v.boxPlot();
       v.treeMap();
       v.heatMap();
       v.statistics();
       v.corrPlot();

    });

    function viewerCheck(viewers:DG.Viewer[], viewer:string): void {
        let check = false;
        for(let i:number = 0; i < viewers.length; i++){
            if (viewers[i].type == viewer){
                check = true;
                break;
            }
        }
        if (check == false){
            throw viewer + 'not found';
        }
    }

    test('Scatter Plot', async () => {
        viewerCheck(Array.from(v.viewers), 'Scatter plot')
    });

    test('BarChart', async () => {
        viewerCheck(Array.from(v.viewers), 'Bar chart')
    });

    test('Histogram', async () => {
        viewerCheck(Array.from(v.viewers), 'Histogram')
    });

    test('PieChart', async () => {
        viewerCheck(Array.from(v.viewers), 'Pie chart')
    });

    test('LineChart', async () => {
        viewerCheck(Array.from(v.viewers), 'Line chart')
    });

    test('TrellisPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'Trellis plot')
    });

    test('MatrixPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'Matrix plot')
    });

    test('ScatterPlot3D', async () => {
        viewerCheck(Array.from(v.viewers), '3d scatter plot')
    });

    test('DensityPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'Density plot')
    });

    test('PCPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'PC Plot')
    });

    test('WordCloud', async () => {
        viewerCheck(Array.from(v.viewers), 'Word cloud')
    });

    test('NetworkDiagram', async () => {
        viewerCheck(Array.from(v.viewers), 'Network diagram')
    });

    test('BoxPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'Box plot')
    });

    test('TreeMap', async () => {
        viewerCheck(Array.from(v.viewers), 'Tree map')
    });

    test('HeatMap', async () => {
        viewerCheck(Array.from(v.viewers), 'Heat map')
    });

    test('Statistics', async () => {
        viewerCheck(Array.from(v.viewers), 'Statistics')
    });

    test('CorrelationPlot', async () => {
        viewerCheck(Array.from(v.viewers), 'Correlation plot')
    });


    after(async () => {
        v.close();
        grok.shell.closeAll();
    });

});