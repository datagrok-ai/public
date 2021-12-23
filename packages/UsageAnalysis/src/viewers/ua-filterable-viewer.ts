import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import {UaFilter} from "../filter2";
import {UaViewer} from "./ua-viewer";
import {UaDataFrameViewer} from "./ua-data-frame-viewer";

export class UaFilterableViewer extends UaDataFrameViewer {

    public constructor(name: string, queryName: string, viewerFunction: Function, setStyle?: Function) {
        super(name, queryName, viewerFunction, setStyle);
    }

    reload(filter: UaFilter) {
        this.root.innerHTML = '';
        this.root.append(
            this.addCardWithFilters(
                this.name,
                this.queryName,
                filter,
                this.viewerFunction
            )
        );
    }

    init(): void {
    }

}