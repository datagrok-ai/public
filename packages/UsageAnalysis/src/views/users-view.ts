import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import {UaQueryViewer} from "../viewers/ua-query-viewer";

export class UsersView extends UaView {

    constructor(uaToolbox: UaToolbox) {
        super('Users', uaToolbox);
    }

    async initViewers(): Promise<void> {
        let uniqueUsersViewer = new UaFilterableViewer(
            'Unique Users',
            'UniqueUsers',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
        );
        this.viewers.push(uniqueUsersViewer);

        let uniqueSessionsViewer = new UaFilterableViewer(
            'Unique Sessions',
            'UniqueSessions',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t, UaQueryViewer.defaultChartOptions).root
        );
        this.viewers.push(uniqueSessionsViewer);

        let usageViewer = new UaFilterableViewer(
            'Usage',
            'Usage',
            (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions).root
        );
        this.viewers.push(usageViewer);

        let topPackagesByUsers = new UaFilterableViewer(
            'Top Packages By Users',
            'TopPackagesByUsers',
            (t: DG.DataFrame) => DG.Viewer.scatterPlot(t, UaQueryViewer.defaultChartOptions).root
        );
        this.viewers.push(topPackagesByUsers);

        this.root.append(ui.divV([
            ui.divH([ui.block([usageViewer.root])]),
            ui.divH([ui.block50([uniqueUsersViewer.root]), ui.block50([uniqueSessionsViewer.root])]),
            ui.divH([ui.block([usageViewer.root])]),
        ]));

    }
}