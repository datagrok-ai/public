import {
    TYPE_COLUMN,
    TYPE_DATA_FRAME, TYPE_EVENT_DATA, TYPE_GRID_CELL_RENDER_ARGS,
    TYPE_MENU,
    TYPE_PROJECT,
    TYPE_PROPERTY,
    TYPE_SEMANTIC_VALUE, TYPE_TABLE_VIEW,
    TYPE_USER, TYPE_VIEW, TYPES_SCALAR
} from "./const";
import {Column, DataFrame} from "./dataframe";
import {Project, Property, SemanticValue, User} from "./entities";
import {Balloon, Menu} from "./ui";
import {GridCellRenderArgs} from "./grid";
import {EventData} from "./events";
import {TableView, View} from "./view";

export function _toJs(d) { return _wrap(d, false); }

/** Instantiates the corresponding JS handler for the Dart object [d]. */
export function _wrap(d, check = true) {
    let type = grok_GetType(d);
    switch (type) {
        case TYPE_DATA_FRAME: return new DataFrame(d);
        case TYPE_COLUMN: return new Column(d);
        case TYPE_PROPERTY: return new Property(d);
        case TYPE_PROJECT: return new Project(d);
        case TYPE_USER: return new User(d);
        case TYPE_SEMANTIC_VALUE: return new SemanticValue(d);
        case TYPE_MENU: return new Menu(d);
        case TYPE_GRID_CELL_RENDER_ARGS: return new GridCellRenderArgs(d);
        case TYPE_EVENT_DATA: return new EventData(d);
        case TYPE_VIEW: return new View(d);
        case TYPE_TABLE_VIEW: return new TableView(d);
    }

    if (check)
        throw `Not supported type: ${type}`;

    return d;
}

export function _toDart(x) {
    return (typeof x.d !== 'undefined') ? x.d : x;
}

export function _toJson(x) {
    return x === null ? null : JSON.stringify(x);
}

export function _jsThen(promise, f) {
    promise.then(f);
}

export function paramsToJs(params) {
    let result = [];
    for (let i = 0; i < params.length; i++) {
        let type = grok_GetType(params[i]);
        if (type !== null && !TYPES_SCALAR.has(type))
            result.push(_wrap(params[i]));
        else
            result.push(params[i]);
    }

    return result;
}

/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/

let time = function(s, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${s}: ${stop - start}ms`);
    Balloon.info(`${s}: ${stop - start}ms`);
    return result;
};
