import {TYPE} from "./const";
import {Column, DataFrame} from "./dataframe";
import {Project, Property, SemanticValue, User} from "./entities";
import * as ui from "../ui";
import {EventData} from "./events";
import {TableView, View} from "./view";
import {Menu} from "./ui_classes";

export let factories = {};

export function _toJs(d) { return _wrap(d, false); }

/** Instantiates the corresponding JS handler for the Dart object [d]. */
export function _wrap(d, check = true) {
    let type = grok_GetType(d);
    switch (type) {
        case TYPE.DATA_FRAME: return new DataFrame(d);
        case TYPE.COLUMN: return new Column(d);
        case TYPE.PROPERTY: return new Property(d);
        case TYPE.PROJECT: return new Project(d);
        case TYPE.USER: return new User(d);
        case TYPE.SEMANTIC_VALUE: return new SemanticValue(d);
        case TYPE.MENU: return new Menu(d);
        case TYPE.GRID_CELL_RENDER_ARGS: return factories["GridCellRenderArgs"](d);
        case TYPE.EVENT_DATA: return new EventData(d);
        case TYPE.VIEW: return new View(d);
        case TYPE.TABLE_VIEW: return new TableView(d);
    }

    if (check)
        throw `Not supported type: ${type}`;

    return d;
}

export function _toDart(x) {
    return (typeof x.d !== 'undefined') ? x.d : x;
}

