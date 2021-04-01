/** A viewer that is typically docked inside a [TableView]. */
import {TYPE, VIEWER, ViewerPropertyType, ViewerType} from "./const";
import {Column, DataFrame} from "./dataframe.js";
import {DateTime, Property} from "./entities";
import {ObjectPropertyBag, Widget} from "./widgets";
import {_toJson} from "./utils";
import {toJs} from "./wrappers";
import {StreamSubscription, __obs} from "./events";
import * as rxjs from "rxjs";
import { Grid } from "./grid";

declare let DG: any;
declare let ui: any;
let api = <any>window;

export class TypedEventArgs {
  d: any;
  constructor(d: any) {
    this.d = d;
  }

  get type(): string {
    return api.grok_TypedEventArgs_Get_Type(this.d);
  }

  get data(): any {
    let data = api.grok_TypedEventArgs_Get_Data(this.d);
    let jsData = toJs(data);
    return jsData;
  }
}

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
export class Viewer extends Widget {
  d: any;
  props: ObjectPropertyBag;

  /** @constructs Viewer */
  constructor(d: any) {
    super();
    this.d = d;

    /** @member {ObjectPropertyBag} */
    this.props = new ObjectPropertyBag(this, api.grok_Viewer_Get_Look(this.d));
  }

  /** Creates a new viewer of the specified type.
   * @param {ViewerType} viewerType
   * @param {DataFrame} table
   * @param options
   * @returns {Viewer} */
  static fromType(viewerType: ViewerType, table: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_FromType(viewerType, table.d, _toJson(options)));
  }

  static getViewerTypes(): ViewerType[] {
    return api.grok_Viewer_GetViewerTypes();
  }

  /** 
   * Sets viewer options. See also {@link getOptions}
   *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot
   *  @param {object} map */
  setOptions(map: {type: ViewerType}): void {
    api.grok_Viewer_Options(this.d, JSON.stringify(map));
  }

  /** 
   * Gets viewer options. See also {@link getOptions}
   *  Sample: https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot
   *  @returns {object} */
  getOptions(): {type: ViewerType} {
    return JSON.parse(api.grok_Viewer_Serialize(this.d));
  }

  getInfo(): object {
    return api.grok_Viewer_GetInfo(this.d);
  }

  getProperties(): Property[] {
    return api.grok_Viewer_Get_Properties(this.d);
  }

  /** Closes and detaches the viewer. */
  close(): void {
    api.grok_Viewer_Close(this.d);
  }

  /** Visual root.
   * @type {HTMLElement} */
  get root(): HTMLElement {
    return api.grok_Viewer_Root(this.d);
  }

  /** Returns viewer type (see VIEWER constants)
   * @returns {string} */
  get type(): ViewerType {
    return this.getOptions().type;
  }

  get table(): DataFrame {
    return toJs(api.grok_Viewer_Get_DataFrame(this.d));
  }

  /** @type {DataFrame} */
  get dataFrame(): DataFrame {
    return toJs(api.grok_Viewer_Get_DataFrame(this.d));
  }

  set dataFrame(t: DataFrame) {
    api.grok_Viewer_Set_DataFrame(this.d, t == null ? null : t.d);
  }

  static grid(t: DataFrame, options: object | null = null): Grid {
    return new DG.Grid(api.grok_Viewer_Grid(t.d, _toJson(options)));
  }

  static histogram(t: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_Histogram(t.d, _toJson(options)));
  }

  static barChart(t: DataFrame, options: object | null = null): Viewer {
    return Viewer.fromType(VIEWER.BAR_CHART, t, options);
  }

  static heatMap(t: DataFrame, options: object | null = null): Viewer {
    return Viewer.fromType(VIEWER.HEAT_MAP, t, options);
  }

  static boxPlot(t: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_BoxPlot(t.d, _toJson(options)));
  }

  static filters(t: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_Filters(t.d, _toJson(options)));
  }

  static scatterPlot(t: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_ScatterPlot(t.d, _toJson(options)));
  }

  static lineChart(t: DataFrame, options: object | null = null): Viewer {
    return new Viewer(api.grok_Viewer_LineChart(t.d, _toJson(options)));
  }

  /** Observes platform events with the specified eventId.
   * @param {string} eventId
   * @returns {Observable} */
  onEvent(eventId: string | null = null): rxjs.Observable<any> {
    if (eventId !== null)
      return __obs(eventId, this.d);

    let dartStream = api.grok_Viewer_Get_EventBus_Events(this.d);
    return rxjs.fromEventPattern(
      function (handler) {
        return api.grok_Stream_Listen(dartStream, function (x: any) {
          handler(new TypedEventArgs(x));
        });
      },
      function (handler, d) {
        new StreamSubscription(d).cancel();
      }
    );
  }
}


/** Subclass JsViewer to implement a DataFrame-bound Datagrok viewer in JavaScript.
 *  See an example on github: {@link https://github.com/datagrok-ai/public/tree/master/packages/Leaflet}
 *  */
export class JsViewer extends Widget {
  dataFrame: DataFrame | null;
  subs: StreamSubscription[];
  obs: rxjs.Observable<any>[];
  props: ObjectPropertyBag;

  /** @constructs JsViewer */
  constructor() {
    super();

    /** @type {HTMLElement} */
    this.root = ui.box();

    /** @type {DataFrame} */
    this.dataFrame = null;

    /** @type {StreamSubscription[]} */
    this.subs = [];  // stream subscriptions - will be canceled when the viewer is detached

    /** @member {Observable[]} */
    this.obs = [];

    /** @member {ObjectPropertyBag} */
    this.props = new ObjectPropertyBag(this);
  }

  onFrameAttached(dataFrame: DataFrame): void {
    this.dataFrame = dataFrame;
    this.onTableAttached();
  }

  /** Gets called when a table is attached to the viewer. */
  onTableAttached(): void {
  }

  /** Gets called when this viewer is detached. */
  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Gets property by name (case-sensitive).
   * @param {string} name
   * @returns {Property} */
  getProperty(name: string): Property | undefined {
    return this.getProperties().find((p) => p.name === name);
  }

  getProperties(): Property[] {
    return this._properties;
  }

  /** cleanup() will get called when the viewer is disposed
   * @param {Function} cleanup */
  registerCleanup(cleanup: Function): void {
    api.grok_Widget_RegisterCleanup(this.root, cleanup);
  }

  /**
   * @param {Observable} observable
   * @returns {Observable} */
  _obs(observable: rxjs.Observable<any>): rxjs.Observable<any> {
    this.obs.push(observable);
    return observable;
  }

  /** Returns the column bound to the specified data property.
   *  Note that "ColumnName" suffix (this determines whether this is a data property) should be omitted.
   * @param {string} dataPropertyName
   * @param {object} options
   * @returns {Column} */
  column(dataPropertyName: string, options: {} | null = null): Column {
    return this.addProperty(`${dataPropertyName}ColumnName`, TYPE.STRING, null, options);
  }

  /** Registers an integer property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {number} defaultValue
   * @param {object} options
   * @returns {number} */
  int(propertyName: ViewerPropertyType, defaultValue: number | null = null, options: {} | null = null): number {
    return this.addProperty(propertyName, TYPE.INT, defaultValue, options);
  }

  /** Registers a floating point property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {number} defaultValue
   * @param {object} options
   * @returns {number} */
  float(propertyName: ViewerPropertyType, defaultValue: number | null = null, options: {} | null = null): number {
    return this.addProperty(propertyName, TYPE.FLOAT, defaultValue, options);
  }

  /** Registers a string property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {string} defaultValue
   * @param {object} options
   * @returns {string} */
  string(propertyName: ViewerPropertyType, defaultValue: string | null = null, options: {} | null = null): string {
    return this.addProperty(propertyName, TYPE.STRING, defaultValue, options);
  }

  /** Registers a string list property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {string[]} defaultValue
   * @param {object} options
   * @returns {string[]} */
  stringList(propertyName: ViewerPropertyType, defaultValue: string[] | null = null, options: {} | null = null): string[] {
    return this.addProperty(propertyName, TYPE.STRING_LIST, defaultValue, options);
  }

  /** Registers a boolean property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {boolean} defaultValue
   * @param {object} options
   * @returns {boolean} */
  bool(propertyName: ViewerPropertyType, defaultValue: boolean | null = null, options: {} | null = null): boolean {
    return this.addProperty(propertyName, TYPE.BOOL, defaultValue, options);
  }

  /** Registers a datetime property with the specified name and defaultValue
   * @param {ViewerPropertyType} propertyName
   * @param {DateTime} defaultValue
   * @param {object} options
   * @returns {DateTime} */
  dateTime(propertyName: ViewerPropertyType, defaultValue: DateTime | null = null, options: {} | null = null): DateTime {
    return this.addProperty(propertyName, TYPE.DATE_TIME, defaultValue, options);
  }
}
