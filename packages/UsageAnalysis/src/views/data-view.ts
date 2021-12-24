import {UaView} from "./ua-view";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class DataView extends UaView {

    constructor(uaToolbox: UaToolbox) {
        super('Data', uaToolbox);
    }

    async initViewers() : Promise<void> {
        let queriesViewer = new UaFilterableViewer(
            'Queries',
            'Queries1',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t).root
        );
        this.viewers.push(queriesViewer);

        let topQueriesViewer = new UaFilterableViewer(
            'Top Queries',
            'TopQueries',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topQueriesViewer);

        let topConnectionsViewer = new UaFilterableViewer(
            'Top Connections',
            'TopConnections',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topConnectionsViewer);

        this.root.append(ui.divV([
            ui.divH([ui.block([queriesViewer.root])]),
            ui.divH([ui.block50([topQueriesViewer.root]), ui.block50([topConnectionsViewer.root])]),
        ]));

    }
}