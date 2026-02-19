/**
 * Base widget classes: TypedEventArgs, ObjectPropertyBag, Widget, DartWidget, DartWrapper.
 * @module widgets/base
 */

import {toDart, toJs} from "../wrappers";
import {Subscription} from "rxjs";
import {Func, Property, IProperty} from "../entities";
import {DataFrame} from "../dataframe";
import {Type} from "../const";
import {MapProxy} from "../proxies";
import {IDartApi} from "../api/grok_api.g";

declare let DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class TypedEventArgs<TData> {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** Event type id */
  get type(): string {
    return api.grok_TypedEventArgs_Get_Type(this.dart);
  }

  get data(): TData {
    let data = api.grok_TypedEventArgs_Get_Data(this.dart);
    return toJs(data);
  }

  /** Event arguments. Only applies when data is EventData */
  get args(): {[key: string]: any} | null {
    // @ts-ignore
    if (!this.data?.dart)
      return null;

    // @ts-ignore
    return api.grok_EventData_Get_Args(this.data.dart);
  }
}


/**
 * A dynamic property bag that provides convenient access to an object's properties.
 *
 * This class uses a JavaScript Proxy to allow accessing properties directly by name:
 * ```typescript
 * viewer.props.showAxes = true;  // Instead of: viewer.setOptions({showAxes: true})
 * console.log(viewer.props.color);
 * ```
 *
 * The Proxy intercepts property access and routes it through the underlying
 * property system, which handles type conversion and validation.
 *
 * **Note**: Due to the dynamic Proxy nature, TypeScript cannot provide
 * compile-time type checking for property access. Use `getProperties()` to
 * discover available properties at runtime.
 *
 * @example
 * // Access viewer properties
 * const scatter = view.scatterPlot();
 * scatter.props.xColumnName = 'age';
 * scatter.props.yColumnName = 'weight';
 *
 * // Enumerate available properties
 * for (const prop of scatter.props.getProperties()) {
 *   console.log(prop.name, prop.propertyType);
 * }
 */
export class ObjectPropertyBag {
  /** The source object whose properties are being wrapped */
  source: any;

  constructor(source: any, x: any = null) {

    /** @member {Object} */
    this.source = source;

    if (x == null)
      x = source;

    let props = this;

    const handler = {
      ownKeys(target: any) {
        return props.getProperties().map((p: Property) => p.name);
      },
      has(target: any, name: string) {
        return props.getProperties().find((p: Property) => p.name === name) !== null;
      },
      getOwnPropertyDescriptor(target: any, key: any) {
        return {
          enumerable: true,
          configurable: true
        };
      },
      set(target: any, name: string, value: any) {
        if (<any>x instanceof DartWidget)
          x = toDart(x);
        target.getProperty(name).set(x, value);
        return true;
      },
      get(target: any, name: any) {
        const own = ['__proto__', 'hasProperty', 'getProperty', 'getProperties',
          'get', 'set', 'setAll', 'apply', 'setDefault', 'resetDefault'];
        if (own.includes(name) || props.hasOwnProperty(name))
          // @ts-ignore
          return props[name];
        if (<any>x instanceof DartWidget)
          x = toDart(x);
        return target.getProperty(name)?.get(x);
      }
    }

    return new Proxy(this, handler);
  }

  /**
   * Gets the value of the specified property
   * @param {string} propertyName
   * @returns {object}
   * */
  get(propertyName: string): object {
    return this.getProperty(propertyName).get(this.source);
  }

  /**
   * Sets the value of the specified property
   * @param {string} propertyName
   * @param {object} propertyValue
   * */
  set(propertyName: string, propertyValue: object) {
    this.getProperty(propertyName).set(this.source, propertyValue);
  }

  /** Sets all properties according to the passed object containing key-value pairs */
  setAll(params: object) {
    for (const [k, v] of Object.entries(params))
      this.set(k, v);
  }

  /** @returns {Property[]} */
  getProperties(): Property[] {
    return this.source.getProperties();
  }

  /** Gets property by name (case-sensitive).
   * @param {string} name
   * @returns {Property} */
  getProperty(name: string): Property {
    let property = this.getProperties().find((p) => p.name === name);
    if (typeof property == 'undefined')
      throw `Property not found: ${name}`;
    return property;
  }

  /**
   * @param {string} name
   * @returns {boolean} */
  hasProperty(name: string): boolean {
    return this.getProperties().findIndex((p) => p.name === name) !== -1;
  }

  /** Sets the current state of viewer properties as the default configuration used to create new viewer
  * instances of this type. Equivalent to the "Pick Up / Apply | Set as Default" context menu command.
  * Read more about viewer commands: https://datagrok.ai/help/visualize/viewers/#common-actions
  * @param data indicates if data settings should be copied.
  * @param style indicates if style (non-data) settings should be copied. */
  setDefault(data: boolean = false, style: boolean = true) {
    if (this.source instanceof DG.Viewer)
      api.grok_Viewer_Props_SetDefault(this.source.dart, data, style);
    else
      throw 'Call failed: object is not Viewer instance';
  }

  static setDefaultProperty(viewerType: string, propertyName: string, propertyValue: any) {
    api.grok_Viewer_Props_SetDefaultProperty(viewerType, propertyName, propertyValue);
  }

  /** Clears the previously remembered default settings for viewers of this type. See also: [setDefault] */
  resetDefault() {
    if (this.source instanceof DG.Viewer)
      api.grok_Viewer_Props_ResetDefault(this.source.dart);
    else
      throw 'Call failed: object is not Viewer instance';
  }
}


/** Base class for controls that have a visual root and a set of properties. */
export class Widget<TSettings = any> {

  public get type(): string { return 'Unknown'; }

  /** Contains auxiliary information */
  public temp: any;

  /** Constructor function. No parameters, returns [Widget]. */
  factory: Func | null = null;

  protected _root: HTMLElement;
  protected _properties: Property[] = [];
  protected _functions: Func[] = [];
  props: TSettings & ObjectPropertyBag; //ObjectPropertyBag;
  subs: Subscription[];
  dart: any;
  isDetached: boolean = false;

  /** @constructs Widget and initializes its root. */
  constructor(widgetRoot: HTMLElement) {
    this._root = widgetRoot;
    // @ts-ignore
    this.props = new ObjectPropertyBag(this);
    this.subs = [];
    this.getProperties = this.getProperties.bind(this);
    this.temp = {};
  }

  /** Returns all currently active widgets. */
  static getAll(): Widget[] {
    return api.grok_Widget_GetAll();
  }

  /** Finds existing widget from its visual root. */
  static find(root: Element): Widget | null {
    return api.grok_Widget_FromRoot(root);
  }

  toDart() {
    if (this.dart == null)
      this.dart = api.grok_Widget_Wrap(this);
    return this.dart;
  }

  /** Registers a subscription to an external event.
   * @param {Subscription} subscription */
  sub(subscription: Subscription): void {
    this.subs.push(subscription);
  }

  /**
   * @param {Object} properties
   * @returns {Widget} */
  apply(properties: object): Widget {
    for (let name of Object.keys(properties))
      if (typeof name !== 'undefined')
        // @ts-ignore
        this.props[name] = properties[name];

     //this.props.apply(properties);
    return this;
  }

  /** Returns all properties of this widget. */
  getProperties(): Property[] { return this._properties; }

  /** Functions that are applicable to this particular widget.
    Used in the UI to display context actions, and for the AI integrations. */
  getFunctions(): Func[] { return this._functions;  }

  /** Gets called when viewer's property is changed.
   * @param {Property} property - or null, if multiple properties were changed. */
  onPropertyChanged(property: Property | null): void {}

  getDartProperties(): any[] {
    return this.getProperties().map((p) => p.dart);
  }

  sourceRowsChanged(): void {};

  onFrameAttached(dataFrame: DataFrame): void {
    if (this.props.hasProperty('dataFrame'))
      this.props.set('dataFrame', dataFrame);
  }

  /** Parent widget up the DOM tree, or null. */
  get parent(): Widget | null { return api.grok_Widget_Get_Parent(this.toDart()); }

  /** Parent widget up the DOM tree, or null. */
  get children(): Widget[] { return api.grok_Widget_Get_Children(this.toDart()); }

  /** Widget's visual root.
   * @type {HTMLElement} */
  get root(): HTMLElement { return this._root; }
  set root(r: HTMLElement) { this._root = r; }

  /** Gets called when a widget is detached and will no longer be used. Typically used for unsubscribing from events.
   * Be sure to call super.detach() if this method is overridden.  */
  detach(): void {
    this.subs.forEach((s) => s.unsubscribe());
    this.isDetached = true;
  }

  /** Registers an property with the specified type, name, and defaultValue.
   *  @see Registered property gets added to {@link properties}.
   *  Returns default value, thus allowing to combine registering a property with the initialization
   *
   * @param {string} propertyName
   * @param {TYPE} propertyType
   * @param defaultValue
   * @param {Object} options
   * @returns {*}
   * @private
   */
  addProperty(propertyName: string, propertyType: Type, defaultValue: any = null, options: { [key: string]: any } & IProperty | null = null): any {
    const fieldName = options?.fieldName ?? propertyName;

    let obj = this;
    // @ts-ignore
    let p = Property.create(propertyName, propertyType, (_) => obj[fieldName], null, defaultValue);
    p.set = function (_: any, x: any) {
      // @ts-ignore
      obj[fieldName] = x;
      obj.onPropertyChanged(p);
    };

    if (options !== null) {
      for (let key of Object.keys(options))
        if (key != 'fieldName')
          api.grok_PropMixin_SetPropertyValue(p.dart, key, options[key]);
    }

    this._properties.push(p);
    return p.defaultValue;
  }

  /** Creates a new widget from the root element. */
  static fromRoot(root: HTMLElement): Widget {
    return new Widget(root);
  }

  // /** Creates a {@see Widget} from the specified React component. */
  // // @ts-ignore
  // static react(reactComponent: React.DOMElement<any, any> | Array<React.DOMElement<any, any>> | React.CElement<any, any> | Array<React.CElement<any, any>> | React.ReactElement | Array<React.ReactElement>): Widget {
  //   let widget = Widget.fromRoot(ui.div());
  //   // @ts-ignore
  //   ReactDOM.render(reactComponent, widget.root);
  //   return widget;
  // }
}


/** JS wrapper for the widget implemented in Dart. */
export class DartWidget extends Widget {
  constructor(dart: any) {
    super(api.grok_Widget_Get_Root(dart));
    this.dart = dart;
    this.temp = new MapProxy(api.grok_Widget_Get_Temp(this.dart));
  }

  get type(): string { return api.grok_Widget_Get_Type(this.dart); }
  get root(): HTMLElement { return api.grok_Widget_Get_Root(this.dart); }
  getProperties(): Property[] { return toJs(api.grok_PropMixin_GetProperties(this.dart)); }
  getFunctions(): Func[] { return toJs(api.grok_Widget_GetFunctions(this.dart)); }
}


export class DartWrapper {
  dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }
}
