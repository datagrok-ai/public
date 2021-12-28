import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {UaFilter} from "../filter2";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableViewer} from "../viewers/ua-filterable-viewer";

export abstract class UaView extends DG.ViewBase {
    uaToolbox: UaToolbox;
    viewers: UaFilterableViewer[] = [];

    protected constructor(viewName: string, uaToolbox: UaToolbox) {
        super();
        this.name = viewName;
        this.uaToolbox = uaToolbox;
        this.toolbox = uaToolbox.rootAccordion.root;

        this.initViewers().then(() => this.initFiltering());
    }

    abstract initViewers() : any; // for sync and async code

    initFiltering() : void {
        this.uaToolbox.filterStream.subscribe((filter) => this.reload(filter));
        this.reload(this.uaToolbox.getFilter())
    }

    reload(filter: UaFilter) : void {
        for (let viewer of this.viewers)
            viewer.reload(filter);
    }
}