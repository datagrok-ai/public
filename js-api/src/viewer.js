/** A viewer that is typically docked inside a [TableView]. */
import {TYPE, VIEWER} from "./const";
import {DataFrame} from "./dataframe.js";
import * as ui from "./../ui.js";
import {Property} from "./entities";
import {Widget} from "./widgets";
import {_onSizeChanged, _toJson} from "./utils";
import {Balloon} from "./widgets";
import {toJs} from "./wrappers";
import {observeStream, StreamSubscription, __obs} from "./events";
import * as rxjs from "rxjs";

export class TypedEventArgs {
  constructor(d) {
    this.d = d;
  }

  get type() {
    return grok_TypedEventArgs_Get_Type(this.d);
  }

  get data() {
    let data = grok_TypedEventArgs_Get_Data(this.d);
    let jsData = toJs(data);
    return jsData;
  }
}

/**
 * Represents a {@link https://datagrok.ai/help/visualize/viewers | viewer}.
 * See also {@link https://datagrok.ai/help/develop/js-api#pre-defined-viewers}
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
export class Viewer {

  /** @constructs Viewer */
  constructor(d) {
    this.d = d;
  }

  /** Creates a new viewer of the specified type.
   * @param {ViewerType} viewerType
   * @param {DataFrame} table
   * @param options
   * @returns {Viewer} */
  static fromType(viewerType, table, options = null) {
    return new Viewer(grok_Viewer_FromType(viewerType, table.d, _toJson(options)));
  }

  static getViewerTypes() {
    return grok_Viewer_GetViewerTypes();
  }

  /** Sets viewer options. See also {@link getOptions}
   *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
   *  @param: {object} map */
  setOptions(map) {
    grok_Viewer_Options(this.d, JSON.stringify(map));
  }

  /** Gets viewer options. See also {@link getOptions}
   *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/scatter-plot
   *  @returns: {object} */
  getOptions() {
    return JSON.parse(grok_Viewer_Serialize(this.d));
  }

  getInfo() {
    return grok_Viewer_GetInfo(this.d);
  }

  get properties() {
    return grok_Viewer_Get_Properties(this.d);
  }

  /** Closes and detaches the viewer. */
  close() {
    grok_Viewer_Close(this.d);
  }

  /** Visual root.
   * @type {HTMLElement} */
  get root() {
    return grok_Viewer_Root(this.d);
  }

  /** Returns viewer type (see VIEWER constants)
   * @returns {string} */
  get type() {
    return this.getOptions().type;
  }

  get table() {
    return toJs(grok_Viewer_Get_DataFrame(this.d));
  }

  /** @type {DataFrame} */
  get dataFrame() {
    return toJs(grok_Viewer_Get_DataFrame(this.d));
  }

  set dataFrame(t) {
    grok_Viewer_Set_DataFrame(this.d, t == null ? null : t.d);
  }

  static grid(t, options = null) {
    return new DG.Grid(grok_Viewer_Grid(t.d, _toJson(options)));
  }

  static histogram(t, options = null) {
    return new Viewer(grok_Viewer_Histogram(t.d, _toJson(options)));
  }

  static barChart(t, options = null) {
    return Viewer.fromType(VIEWER.BAR_CHART, t, options);
  }

  static heatMap(t, options = null) {
    return Viewer.fromType(VIEWER.HEAT_MAP, t, options);
  }

  static boxPlot(t, options = null) {
    return new Viewer(grok_Viewer_BoxPlot(t.d, _toJson(options)));
  }

  static filters(t, options = null) {
    return new Viewer(grok_Viewer_Filters(t.d, _toJson(options)));
  }

  static scatterPlot(t, options = null) {
    return new Viewer(grok_Viewer_ScatterPlot(t.d, _toJson(options)));
  }

  static lineChart(t, options = null) {
    return new Viewer(grok_Viewer_LineChart(t.d, _toJson(options)));
  }

  /** Observes platform events with the specified eventId.
   * @returns {Observable} */
  onEvent(eventId = null) {
    if (eventId !== null)
      return __obs(eventId, this.d);

    let dartStream = grok_Viewer_Get_EventBus_Events(this.d);
    return rxjs.fromEventPattern(
      function (handler) {
        return grok_Stream_Listen(dartStream, function (x) {
          handler(new TypedEventArgs(x));
        });
      },
      function (handler, d) {
        new StreamSubscription(d).cancel();
      }
    );
  }
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

  get(name) {
    return null;
  }
}


export class JsViewerProps {

  constructor(viewer) {

    /** @member {JsViewer} */
    this.viewer = viewer;

    return new Proxy(this, {
      set(target, name, value) {
        this.viewer.getProperty(name).set(null, value);
        return true;
      },
      get(target, name) {
        return this.viewer.getProperty(name).get(null);
      }
    });
  }
}


/** Subclass JsViewer to implement a DataFrame-bound Datagrok viewer in JavaScript.
 *  See an example on github: {@link https://github.com/datagrok-ai/public/tree/master/packages/Leaflet}
 *  */
export class JsViewer extends Widget {

  /** @constructs JsViewer */
  constructor() {
    super();

    /** @type {HTMLElement} */
    this.root = ui.div();

    /** @type {DataFrame} */
    this.dataFrame = null;

    /** @type {StreamSubscription[]} */
    this.subs = [];  // stream subscriptions - will be canceled when the viewer is detached

    /** @member {Observable[]} */
    this.obs = [];

    /** @member {JsViewerProps} */
    this.props = new JsViewerProps(this);
  }

  onFrameAttached(dataFrameHandle) {
    this.dataFrame = new DG.DataFrame(dataFrameHandle);
    this.onTableAttached();
  }

  /** Gets called when a table is attached to the viewer. */
  onTableAttached() {
  }

  /** Gets called when viewer's property is changed.
   * @param {Property} property - or null, if multiple properties were changed. */
  onPropertyChanged(property) {
  }

  /** Gets called when viewer's size is changed.
   * @returns {rxjs.Observable} */
  get onSizeChanged() {
    return _onSizeChanged(this.root);
  }

  /** Gets called when this viewer is detached. */
  detach() {
    //Balloon.info("Detached");
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Gets property by name (case-sensitive).
   * @param {string} name
   * @returns {Property} */
  getProperty(name) {
    return this.properties.find((p) => p.name === name);
  }

  getProperties() {
    return this.properties;
  }

  getDartProperties() {
    return this.getProperties().map((p) => p.d);
  }

  /** cleanup() will get called when the viewer is disposed
   * @param {Function} cleanup */
  registerCleanup(cleanup) {
    grok_Widget_RegisterCleanup(this.root, cleanup);
  }

  /** Registers an property with the specified type, name, and defaultValue.
   *  Registered property gets added to {@see properties}.
   *  Returns default value, thus allowing to combine registering a property with the initialization
   *
   * @param {string} propertyName
   * @param {TYPE} propertyType
   * @param defaultValue
   * @returns {*}
   * @private
   */
  _prop(propertyName, propertyType, defaultValue = null) {
    let obj = this;
    let p = Property.create(propertyName, propertyType, () => obj[propertyName], null, defaultValue);
    p.set = function (_, x) {
      obj[propertyName] = x;
      obj.onPropertyChanged(p);
    };

    this.properties.push(p);
    return p.defaultValue;
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
   * @returns {Column} */
  column(dataPropertyName) {
    return this._prop(`${dataPropertyName}ColumnName`, TYPE.STRING);
  }

  /** Registers an integer property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {number} defaultValue
   * @returns {number} */
  int(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.INT, defaultValue);
  }

  /** Registers a floating point property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {number} defaultValue
   * @returns {number} */
  float(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.FLOAT, defaultValue);
  }

  /** Registers a string property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {string} defaultValue
   * @returns {string} */
  string(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.STRING, defaultValue);
  }

  /** Registers a string list property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {string[]} defaultValue
   * @returns {string[]} */
  stringList(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.STRING_LIST, defaultValue);
  }

  /** Registers a boolean property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {boolean} defaultValue
   * @returns {boolean} */
  bool(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.BOOL, defaultValue);
  }

  /** Registers a datetime property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {DateTime} defaultValue
   * @returns {DateTime} */
  dateTime(propertyName, defaultValue = null) {
    return this._prop(propertyName, TYPE.DATE_TIME, defaultValue);
  }
}
