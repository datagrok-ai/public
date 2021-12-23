import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";

export class EventsView extends UaView {


    constructor(uaToolbox: UaToolbox) {
        super('Events', uaToolbox);
    }

    async initViewers() : Promise<void> {
        let topFunctionsViewer = new UaFilterableViewer(
            'Top Functions',
            'TopFunctions',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topFunctionsViewer);

        let topPackageFunctionsViewer = new UaFilterableViewer(
            'Top Package Functions',
            'TopPackageFunctions',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topPackageFunctionsViewer);

        let topPackagesViewer = new UaFilterableViewer(
            'Top Packages',
            'TopPackages',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topPackagesViewer);

        let eventsViewer = new UaFilterableViewer(
            'Events',
            'Events1',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t).root
        );
        this.viewers.push(eventsViewer);

        let topSourcesViewer = new UaFilterableViewer(
            'Top Sources',
            'TopSources',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        )
        this.viewers.push(topSourcesViewer);

        this.root.append(ui.divV([
                // ui.divH([ui.block([eventsViewer.root])]),
                ui.divH([ui.block50([topPackagesViewer.root]), ui.block50([topPackageFunctionsViewer.root])]),
                ui.divH([ui.block50([topFunctionsViewer.root]), ui.block50([topSourcesViewer.root])])
            ])
        );
    }

}