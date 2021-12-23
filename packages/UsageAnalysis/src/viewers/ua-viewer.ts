import * as ui from "datagrok-api/ui";
import {UaFilter} from "../filter2";

export abstract class UaViewer {
    root: HTMLElement;
    name: string;

    protected constructor(name: string) {
        this.root = ui.div();
        this.name = name;
    }

    abstract init() : void;
    abstract reload(filter: UaFilter) : void;
}
