import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import {TopQueriesUsingDataSource} from "./top-queries-using-data-source";
import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../viewers/ua-query-viewer";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";

export class TopDataSourcesViewer extends UaFilterableViewer {
    public constructor() {
        super(
            'Top Data Sources',
            'TopDataSources',
            (t: DG.DataFrame) => {
                let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
                viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((cats: string[]) => {
                    let viewer = new TopQueriesUsingDataSource(cats[0]);
                    grok.shell.o = ui.block([viewer.root]);
                });
                return viewer.root;
            }
        );
    }

}