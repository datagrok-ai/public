import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";

export class ErrorsView extends UaView {

    constructor(uaToolbox: UaToolbox) {
        super('Errors', uaToolbox);
    }

    async initViewers() : Promise<void> {
        let errorsViewer = new UaFilterableViewer(
            'Errors',
            'Errors1',
            (t: DG.DataFrame) => DG.Viewer.lineChart(t).root
        );
        this.viewers.push(errorsViewer);

        let topErrorsViewer = new UaFilterableViewer(
            'Top Errors',
            'TopErrors',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topErrorsViewer);

        let topErrorSourcesViewer = new UaFilterableViewer(
            'Top Error Sources',
            'TopErrorSources',
            (t: DG.DataFrame) => DG.Viewer.barChart(t, this.defaultBarchartOptions).root
        );
        this.viewers.push(topErrorSourcesViewer);

        this.root.append(ui.divV([
            ui.divH([ui.block([errorsViewer.root])]),
            ui.divH([ui.block50([topErrorsViewer.root]), ui.block50([topErrorSourcesViewer.root])]),
        ]));

    }
}