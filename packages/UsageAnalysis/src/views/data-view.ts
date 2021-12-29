import {UaView} from "./ua-view";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {UaDataFrameViewer} from "../viewers/ua-data-frame-viewer";
import {UaQueryViewer} from "../viewers/ua-query-viewer";
import {TopDataSourcesViewer} from "../drilldown_viewers/top-data-sources-viewer";

export class DataView extends UaView {

    constructor(uaToolbox: UaToolbox) {
        super('Data', uaToolbox);
    }

    async initViewers() : Promise<void> {
        let queriesViewer = new UaFilterableViewer(
            'Queries',
            'Queries1',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
        );
        this.viewers.push(queriesViewer);

        let topQueriesViewer = new UaFilterableViewer(
            'Top Queries',
            'TopQueries',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root
        );
        this.viewers.push(topQueriesViewer);

        let topConnectionsViewer = new UaFilterableViewer(
            'Top Connections',
            'TopConnections',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root
        );
        this.viewers.push(topConnectionsViewer);

        let topDataSourcesViewer = new TopDataSourcesViewer(this.uaToolbox.getFilter());
        this.viewers.push(topDataSourcesViewer);

        this.root.append(ui.divV([
            ui.divH([ui.block([queriesViewer.root])]),
            ui.divH([ui.block50([topQueriesViewer.root]), ui.block50([topConnectionsViewer.root])]),
            ui.divH([ui.block50([topDataSourcesViewer.root])]),
        ]));

    }
}