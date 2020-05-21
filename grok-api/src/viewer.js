/** A viewer that is typically docked inside a [TableView]. */
import {TYPE, VIEWER} from "./const";
import * as ui from "./../ui.js";
import {Property} from "./entities";
import {_toJson} from "./utils";

export class Viewer {
    constructor(d) { this.d = d; }

    static fromType(viewerType, table, options = null) {
        return new Viewer(grok_Viewer_FromType(viewerType, table.d, _toJson(options)));
    }

    options(map) { grok_Viewer_Options(this.d, JSON.stringify(map)); }
    close() { grok_Viewer_Close(this.d); }
    serialize()  { return grok_Viewer_Serialize(this.d); }

    get root() { return grok_Viewer_Root(this.d); }

    static grid        (t, options = null) { return new Viewer(grok_Viewer_Grid(t.d, _toJson(options))); }
    static histogram   (t, options = null) { return new Viewer(grok_Viewer_Histogram(t.d, _toJson(options))); }
    static barChart    (t, options = null) { return Viewer.fromType(VIEWER.BAR_CHART, t, options);  }
    static boxPlot     (t, options = null) { return new Viewer(grok_Viewer_BoxPlot(t.d, _toJson(options)));  }
    static filters     (t, options = null) { return new Viewer(grok_Viewer_Filters(t.d, _toJson(options))); }
    static scatterPlot (t, options = null) { return new Viewer(grok_Viewer_ScatterPlot(t.d, _toJson(options))); }
}

export class JsLookAndFeel {
    constructor() {
        return new Proxy(this, {
            set(target, name, value) {
                target.table.set(name, target.idx, value);
                return true;
            },
            get(target, name) {
                return target.table.get(name, target.idx);
            }
        });
    }

    get(name) { return null; }
}


/** JavaScript implementation of the grok viewer */
export class JsViewer {

    constructor() {
        this.root = ui.div();
        /** @type {Property[]}*/
        this.properties = [];
        this.dataFrameHandle = null;
        this.dataFrame = null;

        this.subs = [];  // stream subscriptions - will be canceled when the viewer is detached

        ui.tools.handleResize(this.root, (w, h) => this.onSizeChanged(w, h));
    }

    onFrameAttached(dataFrameHandle) {}
    onPropertyChanged(property) {}
    onSizeChanged(width, height) {}

    detach() {
        ui.Balloon.info("Detached");
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    /** @returns {Property} */
    getProperty(name) { return this.properties.find((p) => p.name === name); }
    getProperties() { return this.properties; }
    getDartProperties() { return this.getProperties().map((p) => p.d); }

    /** cleanup() will get called when the viewer is disposed **/
    registerCleanup(cleanup) { grok_Widget_RegisterCleanup(this.root, cleanup); }

    _prop(name, type, value = null) {
        let obj = this;
        let p = Property.create(name, type, () => obj[name], null, value);
        p.set = function(_, x) {
            obj[name] = x;
            obj.onPropertyChanged(p);
        };

        this.properties.push(p);
        return p.defaultValue;
    }

    /** @returns {Column} */
    column(name) { return this._prop(`${name}ColumnName`, TYPE.STRING); }

    int(name, value = null) { return this._prop(name, TYPE.INT, value); }
    float(name, value = null) { return this._prop(name, TYPE.FLOAT, value); }
    string(name, value = null) { return this._prop(name, TYPE.STRING, value); }
    bool(name, value = null) { return this._prop(name, TYPE.BOOL, value); }
    dateTime(name, value = null) { return this._prop(name, TYPE.DATE_TIME, value); }
}
