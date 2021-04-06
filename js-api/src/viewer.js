var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
  if (k2 === undefined) k2 = k;
  o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
  Object.defineProperty(o, 'default', { enumerable: true, value: v });
}) : function(o, v) {
  o['default'] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) if (k !== 'default' && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
  __setModuleDefault(result, mod);
  return result;
};
define(['require', 'exports', './const', './widgets', './utils', './wrappers', './events', 'rxjs'], function (require, exports, const_1, widgets_1, utils_1, wrappers_1, events_1, rxjs) {
  'use strict';
  Object.defineProperty(exports, '__esModule', { value: true });
  exports.JsViewer = exports.Viewer = exports.TypedEventArgs = void 0;
  rxjs = __importStar(rxjs);
  let api = window;
  class TypedEventArgs {
    constructor(d) {
      this.d = d;
    }
    get type() {
      return api.grok_TypedEventArgs_Get_Type(this.d);
    }
    get data() {
      let data = api.grok_TypedEventArgs_Get_Data(this.d);
      let jsData = wrappers_1.toJs(data);
      return jsData;
    }
  }
  exports.TypedEventArgs = TypedEventArgs;
  /**
     * Represents a {@link https://datagrok.ai/help/visualize/viewers | viewer}.
     * See also {@link https://datagrok.ai/help/develop/how-to/manipulate-viewers}
     *
     * Use Viewer to control the viewers. To develop a custom viewer, {@see JsViewer}.
     *
     * @example
     * let view = grok.shell.addTableView(grok.data.demo.demog());
     * view.scatterPlot({
         x: 'height',
         y: 'weight',
         size: 'age',
         color: 'race',
       });
     **/
  class Viewer extends widgets_1.Widget {
    /** @constructs Viewer */
    constructor(d) {
      super();
      this.d = d;
      /** @member {ObjectPropertyBag} */
      this.props = new widgets_1.ObjectPropertyBag(this, api.grok_Viewer_Get_Look(this.d));
    }
    /** Creates a new viewer of the specified type.
         * @param {ViewerType} viewerType
         * @param {DataFrame} table
         * @param options
         * @returns {Viewer} */
    static fromType(viewerType, table, options = null) {
      return new Viewer(api.grok_Viewer_FromType(viewerType, table.d, utils_1._toJson(options)));
    }
    static getViewerTypes() {
      return api.grok_Viewer_GetViewerTypes();
    }
    /**
         *  Sets viewer options. See also {@link getOptions}
         *  Sample: {@link https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot}
         *  @param {object} map */
    setOptions(map) {
      api.grok_Viewer_Options(this.d, JSON.stringify(map));
    }
    /**
         * Gets viewer options. See also {@link getOptions}
         *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot
         *  @returns {object} */
    getOptions() {
      return JSON.parse(api.grok_Viewer_Serialize(this.d));
    }
    getInfo() {
      return api.grok_Viewer_GetInfo(this.d);
    }
    getProperties() {
      return api.grok_Viewer_Get_Properties(this.d);
    }
    /** Closes and detaches the viewer. */
    close() {
      api.grok_Viewer_Close(this.d);
    }
    /** Visual root.
         * @type {HTMLElement} */
    get root() {
      return api.grok_Viewer_Root(this.d);
    }
    /** Returns viewer type (see VIEWER constants)
         * @returns {string} */
    get type() {
      return this.getOptions().type;
    }
    get table() {
      return wrappers_1.toJs(api.grok_Viewer_Get_DataFrame(this.d));
    }
    /** @type {DataFrame} */
    get dataFrame() {
      return wrappers_1.toJs(api.grok_Viewer_Get_DataFrame(this.d));
    }
    set dataFrame(t) {
      api.grok_Viewer_Set_DataFrame(this.d, t == null ? null : t.d);
    }
    static grid(t, options = null) {
      return new DG.Grid(api.grok_Viewer_Grid(t.d, utils_1._toJson(options)));
    }
    static histogram(t, options = null) {
      return new Viewer(api.grok_Viewer_Histogram(t.d, utils_1._toJson(options)));
    }
    static barChart(t, options = null) {
      return Viewer.fromType(const_1.VIEWER.BAR_CHART, t, options);
    }
    static heatMap(t, options = null) {
      return Viewer.fromType(const_1.VIEWER.HEAT_MAP, t, options);
    }
    static boxPlot(t, options = null) {
      return new Viewer(api.grok_Viewer_BoxPlot(t.d, utils_1._toJson(options)));
    }
    static filters(t, options = null) {
      return new Viewer(api.grok_Viewer_Filters(t.d, utils_1._toJson(options)));
    }
    static scatterPlot(t, options = null) {
      return new Viewer(api.grok_Viewer_ScatterPlot(t.d, utils_1._toJson(options)));
    }
    static lineChart(t, options = null) {
      return new Viewer(api.grok_Viewer_LineChart(t.d, utils_1._toJson(options)));
    }
    /** Observes platform events with the specified eventId.
         * @param {string} eventId
         * @returns {Observable} */
    onEvent(eventId = null) {
      if (eventId !== null)
        return events_1.__obs(eventId, this.d);
      let dartStream = api.grok_Viewer_Get_EventBus_Events(this.d);
      return rxjs.fromEventPattern(function (handler) {
        return api.grok_Stream_Listen(dartStream, function (x) {
          handler(new TypedEventArgs(x));
        });
      }, function (handler, d) {
        new events_1.StreamSubscription(d).cancel();
      });
    }
  }
  exports.Viewer = Viewer;
  /** Subclass JsViewer to implement a DataFrame-bound Datagrok viewer in JavaScript.
     *  See an example on github: {@link https://github.com/datagrok-ai/public/tree/master/packages/Leaflet}
     *  */
  class JsViewer extends widgets_1.Widget {
    /** @constructs JsViewer */
    constructor() {
      super();
      /** @type {HTMLElement} */
      this.root = ui.box();
      /** @type {DataFrame} */
      this.dataFrame = null;
      /** @type {StreamSubscription[]} */
      this.subs = []; // stream subscriptions - will be canceled when the viewer is detached
      /** @member {Observable[]} */
      this.obs = [];
      /** @member {ObjectPropertyBag} */
      this.props = new widgets_1.ObjectPropertyBag(this);
    }
    onFrameAttached(dataFrame) {
      this.dataFrame = dataFrame;
      this.onTableAttached();
    }
    /** Gets called when a table is attached to the viewer. */
    onTableAttached() {
    }
    /** Gets called when this viewer is detached. */
    detach() {
      this.subs.forEach((sub) => sub.unsubscribe());
    }
    /** Gets property by name (case-sensitive).
         * @param {string} name
         * @returns {Property} */
    getProperty(name) {
      return this.getProperties().find((p) => p.name === name);
    }
    getProperties() {
      return this._properties;
    }
    /** cleanup() will get called when the viewer is disposed
         * @param {Function} cleanup */
    registerCleanup(cleanup) {
      api.grok_Widget_RegisterCleanup(this.root, cleanup);
    }
    /**
         * @param {Observable} observable
         * @returns {Observable} */
    _obs(observable) {
      this.obs.push(observable);
      return observable;
    }
    /** Returns the column bound to the specified data property.
         *  Note that "ColumnName" suffix (this determines whether this is a data property) should be omitted.
         * @param {string} dataPropertyName
         * @param {object} options
         * @returns {Column} */
    column(dataPropertyName, options = null) {
      return this.addProperty(`${dataPropertyName}ColumnName`, const_1.TYPE.STRING, null, options);
    }
    /** Registers an integer property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {number} defaultValue
         * @param {object} options
         * @returns {number} */
    int(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.INT, defaultValue, options);
    }
    /** Registers a floating point property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {number} defaultValue
         * @param {object} options
         * @returns {number} */
    float(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.FLOAT, defaultValue, options);
    }
    /** Registers a string property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {string} defaultValue
         * @param {object} options
         * @returns {string} */
    string(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.STRING, defaultValue, options);
    }
    /** Registers a string list property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {string[]} defaultValue
         * @param {object} options
         * @returns {string[]} */
    stringList(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.STRING_LIST, defaultValue, options);
    }
    /** Registers a boolean property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {boolean} defaultValue
         * @param {object} options
         * @returns {boolean} */
    bool(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.BOOL, defaultValue, options);
    }
    /** Registers a datetime property with the specified name and defaultValue
         * @param {ViewerPropertyType} propertyName
         * @param {DateTime} defaultValue
         * @param {object} options
         * @returns {DateTime} */
    dateTime(propertyName, defaultValue = null, options = null) {
      return this.addProperty(propertyName, const_1.TYPE.DATE_TIME, defaultValue, options);
    }
  }
  exports.JsViewer = JsViewer;
});
