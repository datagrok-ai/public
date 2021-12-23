import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../../css/usage_analysis.css';
import {UaToolbox} from "../ua-toolbox";
import {UaView} from "./ua-view";
import {UaFilter} from "../filter2";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";
import {UaDataFrameViewer} from "../viewers/ua-data-frame-viewer";
import {UaViewer} from "../viewers/ua-viewer";

export class OverviewView extends UaView {

    constructor(uaToolbox: UaToolbox) {
        super('Overview', uaToolbox);
    }

    async initViewers() : Promise<void> {
        let servicesViewer = await this.getServices();

        let uniqueUsersViewer = new UaFilterableViewer(
            'Unique Users',
            'UniqueUsers',
            (t: DG.DataFrame) => {
                let viewer = DG.Viewer.lineChart(t, UaFilterableViewer.splineStyle).root;
                viewer.style.maxHeight = '150px';
                return viewer;
            }
        );
        this.viewers.push(uniqueUsersViewer);

        let eventsViewer = new UaFilterableViewer(
            'Events',
            'Events1',
            (t: DG.DataFrame) => {
                let viewer = DG.Viewer.lineChart(t, UaFilterableViewer.splineStyle).root;
                viewer.style.maxHeight = '150px';
                return viewer;
            }
        );
        this.viewers.push(eventsViewer);

        let errorsViewer = new UaFilterableViewer(
            'Errors',
            'Errors1',
            (t: DG.DataFrame) => {
                let viewer = DG.Viewer.lineChart(t, UaFilterableViewer.splineStyle).root;
                viewer.style.maxHeight = '150px';
                return viewer;
            }
        );
        this.viewers.push(errorsViewer);

        let uniqueUsersListViewer = new UaFilterableViewer(
            'Unique Users List',
            'UniqueUsersList',
            (t: DG.DataFrame) => {
                let ids = Array.from(t.getCol('id').values());
                return ui.wait(async () =>  ui.list(await grok.dapi.getEntities(ids)));
            },
            (host: HTMLElement) => {
                host.style.overflow="auto";
                host.style.height="94.5%";
            }
        );
        this.viewers.push(uniqueUsersListViewer);

        let totalUsersViewer = new UaDataFrameViewer(
            'Total Users',
            'TotalUsersAndGroups',
            (t: DG.DataFrame) => {
                let list = [
                    ['Total users', t.get('user_count', 0)],
                    ['Total groups', t.get('group_count', 0)]
                ];

                return ui.div([ui.table(list, (item, idx) =>
                    [`${item[0]}:`, item[1]]
                )]);
            }
        );
        this.viewers.push(totalUsersViewer);

        this.root.append(ui.block25([
            ui.block([totalUsersViewer.root]),
            ui.block([servicesViewer]),
            ui.block([uniqueUsersListViewer.root]),
        ]));

        this.root.append(ui.block75([
            ui.block([uniqueUsersViewer.root]),
            ui.block([eventsViewer.root]),
            ui.block([errorsViewer.root])
        ]));
    }


    async getServices() : Promise<HTMLElement> {
        let root = ui.div();
        root.appendChild(ui.h2('Services'));
        let serviceInfos = await grok.dapi.admin.getServiceInfos();
        root.appendChild(ui.table(serviceInfos, (item, idx) =>
            //@ts-ignore
            [`${item.key}:`, $(ui.span([item.status])).addClass(`grok-plugin-status-${item.status.toLowerCase()}`)[0]]));
        return root;
    }

}