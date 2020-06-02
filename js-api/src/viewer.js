/** A viewer that is typically docked inside a [TableView]. */
import {TYPE, VIEWER} from "./const";
import * as ui from "./../ui.js";
import {Property} from "./entities";
import {_toJson} from "./utils";
import {Balloon} from "./ui_classes";

/**
 * Represents a {@link https://datagrok.ai/help/visualize/viewers | viewer}.
 * See also {@link https://datagrok.ai/help/develop/js-api#pre-defined-viewers}
 *
 * Use Viewer to control the viewers. To develop a custom viewer, {@see JsViewer}.
 *
 * @example
 * let view = grok.shell.addTableView(grok.data.testData('demog', 5000));
 * view.scatterPlot({
     x: 'height',
     y: 'weight',
     size: 'age',
     color: 'race',
   });
 **/
export class Viewer {

    /** @constructs Viewer*/
    constructor(d) { this.d = d; }

    /** Creates a new viewer of the specified type.
     * @param {ViewerType} viewerType
     * @param {DataFrame} table
     * @param options
     * @returns {Viewer} */
    static fromType(viewerType, table, options = null) {
        return new Viewer(grok_Viewer_FromType(viewerType, table.d, _toJson(options)));
    }

    /** Sets viewer options. See also {@link getOptions}
     *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
     *  @param: {object} map */
    setOptions(map) { grok_Viewer_Options(this.d, JSON.stringify(map)); }

    /** Gets viewer options. See also {@link getOptions}
     *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
     *  @returns: {object} */
    getOptions()  { return grok_Viewer_Serialize(this.d); }

    /** Closes and detaches the viewer. */
    close() { grok_Viewer_Close(this.d); }

    /** Visual root.
     * @type {HtmlElement} */
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


/** Subclass JsViewer to implement a DataFrame-bound Datagrok viewer in JavaScript.
 *  See an example on github: {@link https://github.com/datagrok-ai/public/tree/master/packages/Leaflet}
 *  */
export class JsViewer {

    /** @constructs JsViewer */
    constructor() {
        this.root = ui.div();

        /** @type {Property[]}*/
        this.properties = [];

        this.dataFrameHandle = null;

        /** @type {DataFrame} */
        this.dataFrame = null;

        /** @type {StreamSubscription[]} */
        this.subs = [];  // stream subscriptions - will be canceled when the viewer is detached

        ui.tools.handleResize(this.root, (w, h) => this.onSizeChanged(w, h));
    }

    /** Gets called when a table is attached to the viewer.
     * @param {DataFrame} dataFrameHandle;
     * */
    onFrameAttached(dataFrameHandle) {}

    onPropertyChanged(property) {}

    onSizeChanged(width, height) {}

    detach() {
        Balloon.info("Detached");
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    /** Gets property ba name (case-sensitive).
     * @param {string} name
     * @returns {Property} */
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
