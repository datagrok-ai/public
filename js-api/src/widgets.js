import { toDart, toJs } from "./wrappers";
import { __obs, _sub, observeStream } from "./events";
import { Property } from "./entities";
let api = window;
export class ObjectPropertyBag {
  constructor(source, x = null) {
    /** @member {Object} */
    this.source = source;
    if (x == null)
      x = source;
    let props = this;
    const handler = {
      ownKeys(target) {
        return props.getProperties().map((p) => p.name);
      },
      has(target, name) {
        return props.getProperties().find((p) => p.name === name) !== null;
      },
      getOwnPropertyDescriptor(target, key) {
        return {
          enumerable: true,
          configurable: true
        };
      },
      set(target, name, value) {
        target.getProperty(name).set(x, value);
        return true;
      },
      get(target, name) {
        const own = ['__proto__', 'hasProperty', 'getProperty', 'getProperties', 'get', 'set', 'apply'];
        if (own.includes(name) || props.hasOwnProperty(name))
          // @ts-ignore
          return props[name];
        return target.getProperty(name).get(x);
      }
    };
    return new Proxy(this, handler);
  }
  // /**
  //  * @param {Object} properties  */
  // apply(properties) {
  //   for (let name of Object.keys(properties)) {
  //     if (this.hasProperty(name))
  //       this.set(name, properties.name);
  //   }
  // }
  /**
   * Gets the value of the specified property
   * @param {string} propertyName
   * @returns {object}
   * */
  get(propertyName) {
    return this.getProperty(propertyName).get(this.source);
  }
  /**
   * Sets the value of the specified property
   * @param {string} propertyName
   * @param {object} propertyValue
   * */
  set(propertyName, propertyValue) {
    this.getProperty(propertyName).set(this.source, propertyValue);
  }
  /** @returns {Property[]} */
  getProperties() {
    return this.source.getProperties();
  }
  /** Gets property by name (case-sensitive).
   * @param {string} name
   * @returns {Property} */
  getProperty(name) {
    console.log(`get ${name}`);
    let property = this.getProperties().find((p) => p.name === name);
    if (typeof property == 'undefined')
      throw `Property not found: ${name}`;
    return property;
  }
  /**
   * @param {string} name
   * @returns {boolean} */
  hasProperty(name) {
    return this.getProperties().findIndex((p) => p.name === name) !== -1;
  }
}
/** Base class for controls that have a visual root and a set of properties. */
export class Widget {
  /** @constructs Widget and initializes its root. */
  constructor(widgetRoot = null) {
    /** @member {HTMLElement} */
    this._root = widgetRoot;
    /** @member {Property[]}*/
    this._properties = [];
    /** @member {ObjectPropertyBag} */
    this.props = new ObjectPropertyBag(this);
    /** @member {StreamSubscription[]} */
    this.subs = [];
    this.getProperties = this.getProperties.bind(this);
  }
  /** Registers a subscription to an external event.
   * @param {StreamSubscription} subscription */
  sub(subscription) {
    this.subs.push(subscription);
  }
  /**
   * @param {Object} properties
   * @returns {Widget} */
  apply(properties) {
    for (let name of Object.keys(properties))
      if (typeof name !== 'undefined')
        // @ts-ignore
        this.props[name] = properties[name];
    //this.props.apply(properties);
    return this;
  }
  getProperties() { return this._properties; }
  /** Gets called when viewer's property is changed.
   * @param {Property} property - or null, if multiple properties were changed. */
  onPropertyChanged(property) { }
  getDartProperties() {
    return this.getProperties().map((p) => p.d);
  }
  onFrameAttached(dataFrame) {
    if (this.props.hasProperty('dataFrame'))
      this.props.set('dataFrame', dataFrame);
  }
  /** Widget's visual root.
   * @type {HTMLElement} */
  get root() { return this._root; }
  set root(r) { this._root = r; }
  /** Gets called when a widget is detached and will no longer be used. Typically used for unsubscribing from events.
   * Be sure to call super.detach() if this method is overridden.  */
  detach() {
    this.subs.forEach((s) => s.unsubscribe());
  }
  /** Registers an property with the specified type, name, and defaultValue.
   *  Registered property gets added to {@see properties}.
   *  Returns default value, thus allowing to combine registering a property with the initialization
   *
   * @param {string} propertyName
   * @param {TYPE} propertyType
   * @param defaultValue
   * @param {Object} options
   * @returns {*}
   * @private
   */
  addProperty(propertyName, propertyType, defaultValue = null, options = null) {
    let obj = this;
    // @ts-ignore
    let p = Property.create(propertyName, propertyType, () => obj[propertyName], null, defaultValue);
    p.set = function (_, x) {
      // @ts-ignore
      obj[propertyName] = x;
      obj.onPropertyChanged(p);
    };
    if (options !== null) {
      for (let key of Object.keys(options))
        api.grok_PropMixin_SetPropertyValue(p.d, key, options[key]);
    }
    this._properties.push(p);
    return p.defaultValue;
  }
  /** @returns {Widget} */
  static fromRoot(root) {
    let w = new Widget();
    w.root = root;
    return w;
  }
  /** Creates a {@see Widget} from the specified React component. */
  // @ts-ignore
  static react(reactComponent) {
    let widget = Widget.fromRoot(ui.div());
    // @ts-ignore
    ReactDOM.render(reactComponent, widget.root);
    return widget;
  }
}
/** Base class for DataFrame-bound filtering controls */
export class Filter extends Widget {
  constructor() {
    super(ui.div());
    /** @member {DataFrame} */
    this.dataFrame = null;
  }
  /** Gets called when a data frame is attached.
   * @param {DataFrame} dataFrame*/
  attach(dataFrame) { }
  detach() { }
}
export class DartWidget extends Widget {
  constructor(d) {
    super();
    this.d = d;
  }
  get root() {
    return api.grok_Widget_Get_Root(this.d);
  }
}
/**
 * Accordion control with collapsible/expandable panes.
 * Samples: {@link https://public.datagrok.ai/js/samples/ui/accordion}
 * @extends {DartWidget}
 * */
export class Accordion extends DartWidget {
  /** @constructs Accordion */
  constructor(d) {
    super(d);
  }
  /** Creates a new instance of Accordion */
  static create(key = null) {
    return toJs(api.grok_Accordion(key));
  }
  /** @type {AccordionPane[]} */
  get panes() {
    return api.grok_TabControlBase_Get_Panes(this.d).map(toJs);
  }
  /** Returns a pane with the specified name.
   * @param {string} name
   * @returns {AccordionPane} */
  getPane(name) {
    return toJs(api.grok_TabControlBase_GetPane(this.d, name));
  }
  addPane(name, getContent, expanded = false, before = null) {
    return toJs(api.grok_Accordion_AddPane(this.d, name, getContent, expanded, before !== null ? before.d : null));
  }
}
/** A pane in the {@link Accordion} control. */
export class AccordionPane {
  constructor(d) {
    this.d = d;
  }
  /** Expanded state
   * @type {boolean} */
  get expanded() {
    return api.grok_AccordionPane_Get_Expanded(this.d);
  }
  set expanded(v) {
    api.grok_AccordionPane_Set_Expanded(this.d, v);
  }
  /** @type {string} */
  get name() {
    return api.grok_AccordionPane_Get_Name(this.d);
  }
  set name(name) {
    api.grok_AccordionPane_Set_Name(this.d, name);
  }
}
export class TabControl {
  constructor(d) {
    this.d = d;
  }
  static create(vertical = false) {
    return toJs(api.grok_TabControl(vertical));
  }
  get root() {
    return api.grok_Widget_Get_Root(this.d);
  }
  get header() {
    return api.grok_TabControlBase_Get_Header(this.d);
  }
  get panes() {
    return api.grok_TabControlBase_Get_Panes(this.d).map(toJs);
  }
  getPane(name) {
    return toJs(api.grok_TabControlBase_GetPane(this.d, name));
  }
  addPane(name, getContent, icon = null) {
    return toJs(api.grok_TabControlBase_AddPane(this.d, name, getContent, icon));
  }
  clear() {
    api.grok_TabControlBase_Clear(this.d);
  }
  get currentPane() {
    return api.grok_TabControlBase_Get_CurrentPane(this.d);
  }
  set currentPane(v) {
    api.grok_TabControlBase_Set_CurrentPane(this.d, v.d);
  }
}
export class TabPane {
  constructor(d) {
    this.d = d;
  }
  get expanded() {
    return api.grok_AccordionPane_Get_Expanded(this.d);
  }
  set expanded(v) {
    api.grok_AccordionPane_Set_Expanded(this.d, v);
  }
  get name() {
    return api.grok_AccordionPane_Get_Name(this.d);
  }
  set name(name) {
    api.grok_AccordionPane_Set_Name(this.d, name);
  }
}
export class ToolboxPage {
  constructor(d) {
    this.d = d;
  }
  get accordion() {
    return toJs(api.grok_ToolboxPage_Get_Accordion(this.d));
  }
}
/**
 * A non-modal dialog.
 * Sample: https://public.datagrok.ai/js/samples/ui/dialogs
 *
 * @example
 * ui.dialog('Windows')
 *   .add(ui.)
 *   .add(ui.span(['People of Earth, your attention, please… ']))
 *   .onOK(() => { grok.shell.info('OK!'); })
 *   .show();
 * */
export class Dialog {
  constructor(d) {
    this.d = d;
  }
  static create(title = '') {
    return new Dialog(api.grok_Dialog(title));
  }
  /**
   * Sets the OK button handler, and shows the OK button
   * @param {Function} handler
   * @returns {Dialog} */
  onOK(handler) {
    api.grok_Dialog_OnOK(this.d, handler);
    return this;
  }
  /**
   * Sets the CANCEL button handler
   * @param {Function} handler
   * @returns {Dialog} */
  onCancel(handler) {
    api.grok_Dialog_OnCancel(this.d, handler);
    return this;
  }
  /** @returns {Observable} */
  get onClose() {
    return __obs('d4-dialog-closed', this.d);
  }
  // Using __obs is a recommended method. The below are obsolete and shall not be used:
  // onClose(handler) { api.grok_Dialog_OnClose(this.d, handler); return this; }
  // onClose(handler) { let s = _sub(api.grok_Dialog_OnClose(this.d, () => { handler(); s.cancel(); })); return this; }
  /** @returns {Dialog}
   * @param {{modal: boolean, fullScreen: boolean, center: boolean, centerAt: Element, x: number, y: number, width: number, height: number}|{}} options
   * */
  show(options) {
    options = options || {};
    api.grok_Dialog_Show(this.d, options.modal, options.fullScreen, options.center, options.centerAt, options.x, options.y, options.width, options.height);
    return this;
  }
  /** @returns {Dialog}
   * @param {boolean} fullScreen  */
  showModal(fullScreen) {
    api.grok_Dialog_Show(this.d, true, fullScreen, false, null, null, null, null, null);
    return this;
  }
  /** Adds content to the dialog.
   * @param {HTMLElement | Widget | InputBase} content
   * @returns {Dialog} */
  add(content) {
    api.grok_Dialog_Add(this.d, toDart(ui.extract(content)));
    return this;
  }
  /** Closes the dialog. */
  close() {
    api.grok_Dialog_Close(this.d);
  }
  /** Returns command button with the specified text.
   * @param {string} text
   * @returns {HTMLButtonElement}
   * */
  getButton(text) {
    return api.grok_Dialog_GetButton(this.d, text);
  }
  /** Adds command button with the specified text.
   * @param {string} text
   * @param {Function} action
   * @returns {Dialog}
   * */
  addButton(text, action, index = 0, tooltip = null) {
    api.grok_Dialog_AddButton(this.d, text, action, index, tooltip);
    return this;
  }
  /** Adds context action with the specified text.
   * @param {string} text
   * @param {Function} action
   * @returns {Dialog}
   * */
  addContextAction(text, action) {
    api.grok_Dialog_AddContextAction(this.d, text, action);
    return this;
  }
}
/**
 * Menu (either top menu or popup menu).
 * Top menu sample: {@link https://public.datagrok.ai/js/samples/ui/menu}
 * Popup menu sample: {@link https://public.datagrok.ai/js/samples/ui/popup-menu}
 *
 * @example
 * DG.Menu.popup()
 *   .item('Show info', () => grok.shell.info('Info'))
 *   .separator()
 *   .items(['First', 'Second'], showBalloon)
 *   .show();
 * */
export class Menu {
  constructor(d) {
    this.d = d;
  }
  static create() {
    return toJs(api.grok_Menu());
  }
  /** Creates a popup menu.
   * @returns {Menu} */
  static popup() {
    return toJs(api.grok_Menu_Context());
  }
  /** Finds a child menu item with the specified text.
   * @param {string} text
   * @returns {Menu} */
  find(text) {
    return toJs(api.grok_Menu_Find(this.d, text));
  }
  /** Removes a child menu item with the specified text.
   * @param {string} text */
  remove(text) {
    api.grok_Menu_Remove(this.d, text);
  }
  /** Removes all child menu items. */
  clear() {
    api.grok_Menu_Clear(this.d);
  }
  /** Adds a menu group with the specified text.
   * @param {string} text
   * @returns {Menu} */
  group(text) {
    return toJs(api.grok_Menu_Group(this.d, text));
  }
  /** Adds a menu group with the specified text and handler.
   * @param {string} text
   * @param {Function} onClick - callback with no parameters
   * @returns {Menu} */
  item(text, onClick) {
    return toJs(api.grok_Menu_Item(this.d, text, onClick));
  }
  /** For each item in items, adds a menu group with the specified text and handler.
   * Returns this.
   * @param {string[]} items
   * @param {Function} onClick - a callback with one parameter
   * @returns {Menu} */
  items(items, onClick) {
    return toJs(api.grok_Menu_Items(this.d, items, onClick));
  }
  /** Adds a separator line.
   *  @returns {Menu} */
  separator() {
    return toJs(api.grok_Menu_Separator(this.d));
  }
  /** Shows the menu.
   * @returns {Menu} */
  show() {
    return toJs(api.grok_Menu_Show(this.d));
  }
}
/** Balloon-style visual notifications. */
export class Balloon {
  /** Shows information message (green background) */
  info(s) {
    api.grok_Balloon(s, 'info');
  }
  /** Shows information message (red background) */
  error(s) {
    api.grok_Balloon(s, 'error');
  }
}
/** Input control base. Could be used for editing {@link Property} values as well. */
export class InputBase {
  constructor(d, onChanged = null) {
    this.d = d;
    if (onChanged != null)
      this.onChanged((_) => onChanged(this.value));
  }
  get root() {
    return api.grok_InputBase_Get_Root(this.d);
  }
  ;
  get caption() {
    return api.grok_InputBase_Get_Caption(this.d);
  }
  ;
  get format() {
    return api.grok_InputBase_Get_Format(this.d);
  }
  ;
  get captionLabel() {
    return api.grok_InputBase_Get_CaptionLabel(this.d);
  }
  ;
  get input() {
    return api.grok_InputBase_Get_Input(this.d);
  }
  ;
  get nullable() {
    return api.grok_InputBase_Get_Nullable(this.d);
  }
  ;
  set nullable(v) {
    api.grok_InputBase_Set_Nullable(this.d, v);
  }
  ;
  get value() {
    return toJs(api.grok_InputBase_Get_Value(this.d));
  }
  ;
  set value(x) {
    toDart(api.grok_InputBase_Set_Value(this.d, x));
  }
  ;
  get stringValue() {
    return api.grok_InputBase_Get_StringValue(this.d);
  }
  ;
  set stringValue(s) {
    api.grok_InputBase_Set_StringValue(this.d, s);
  }
  ;
  get readOnly() {
    return api.grok_InputBase_Get_ReadOnly(this.d);
  }
  ;
  set readOnly(v) {
    api.grok_InputBase_Set_ReadOnly(this.d, v);
  }
  ;
  get enabled() {
    return api.grok_InputBase_Get_Enabled(this.d);
  }
  ;
  set enabled(v) {
    api.grok_InputBase_Set_Enabled(this.d, v);
  }
  ;
  /// Occurs when [value] is changed, either by user or programmatically.
  onChanged(callback) {
    return _sub(api.grok_InputBase_OnChanged(this.d, callback));
  }
  /// Occurs when [value] is changed by user.
  onInput(callback) {
    return _sub(api.grok_InputBase_OnInput(this.d, callback));
  }
  save() {
    return api.grok_InputBase_Save(this.d);
  }
  ;
  load(s) {
    return api.grok_InputBase_Load(this.d, s);
  }
  ;
  init() {
    return api.grok_InputBase_Init(this.d);
  }
  ;
  fireChanged() {
    return api.grok_InputBase_FireChanged(this.d);
  }
  ;
  addCaption(caption) {
    api.grok_InputBase_AddCaption(this.d, caption);
  }
  ;
  addPatternMenu(pattern) {
    api.grok_InputBase_AddPatternMenu(this.d, pattern);
  }
  setTooltip(msg) {
    api.grok_InputBase_SetTooltip(this.d, msg);
  }
  ;
  static forProperty(property) {
    return toJs(api.grok_InputBase_ForProperty(property.d));
  }
}
export class ProgressIndicator {
  constructor(d) {
    this.d = d;
  }
  static create() {
    return toJs(api.grok_ProgressIndicator_Create());
  }
  get percent() {
    return api.grok_ProgressIndicator_Get_Percent(this.d);
  }
  get description() {
    return api.grok_ProgressIndicator_Get_Description(this.d);
  }
  set description(s) {
    api.grok_ProgressIndicator_Set_Description(this.d, s);
  }
  update(percent, description) {
    api.grok_ProgressIndicator_Update(this.d, percent, description);
  }
  log(line) {
    api.grok_ProgressIndicator_Log(this.d, line);
  }
  get onProgressUpdated() {
    return observeStream(api.grok_Progress_Updated(this.d));
  }
  get onLogUpdated() {
    return observeStream(api.grok_Progress_Log_Updated(this.d));
  }
}
export class TaskBarProgressIndicator extends ProgressIndicator {
  static create(name) {
    return toJs(api.grok_TaskBarProgressIndicator_Create(name));
  }
  close() {
    return api.grok_TaskBarProgressIndicator_Close(this.d);
  }
}
export class TagEditor {
  constructor(d) {
    this.d = d;
  }
  static create() {
    return toJs(api.grok_TagEditor());
  }
  get root() {
    return api.grok_TagEditor_Get_Root(this.d);
  }
  get tags() {
    return api.grok_TagEditor_Get_Tags(this.d);
  }
  addTag(tag, notify = true) {
    return api.grok_TagEditor_AddTag(this.d, tag, notify);
  }
  removeTag(tag) {
    api.grok_TagEditor_RemoveTag(this.d, tag);
  }
  clearTags() {
    api.grok_TagEditor_ClearTags(this.d);
  }
  set acceptsDragDrop(predicate) {
    // @ts-ignore
    api.grok_TagEditor_Set_AcceptsDragDrop(this.d, (x) => predicate(toJs(x, false)));
  }
  ;
  set doDrop(action) {
    // @ts-ignore
    api.grok_TagEditor_Set_DoDrop(this.d, (x) => action(toJs(x, false)));
  }
  onChanged(callback) {
    return _sub(api.grok_TagEditor_OnChanged(this.d, callback));
  }
}
export class TagElement {
  constructor(d) {
    this.d = d;
  }
  get tag() {
    return api.grok_TagElement_Get_Tag(this.d);
  }
  ;
  set tag(x) {
    api.grok_TagElement_Set_Tag(this.d, x);
  }
  ;
}
/** Color-related routines. */
export class Color {
  static r(c) {
    return (c >> 16) & 0xFF;
  }
  static g(c) {
    return (c >> 8) & 0xFF;
  }
  static b(c) {
    return c & 0xFF;
  }
  /** Returns i-th categorical color (looping over the palette if needed) */
  static getCategoricalColor(i) {
    return Color.categoricalPalette[i % Color.categoricalPalette.length];
  }
  /** Returns either black or white color, depending on which one would be most contrast to the specified [color] */
  static getContrastColor(color) {
    return api.grok_Color_GetContrastColor(color);
  }
  static toRgb(color) {
    return color === null ? '' : `rgb(${Color.r(color)},${Color.g(color)},${Color.b(color)})`;
  }
  static get categoricalPalette() {
    return api.grok_Color_CategoricalPalette();
  }
  static scale(x, min, max) {
    return min === max ? min : (x - min) / (max - min);
  }
  static get gray() {
    return 0xFF808080;
  }
  static get lightLightGray() {
    return 0xFFF0F0F0;
  }
  static get lightGray() {
    return 0xFFD3D3D3;
  }
  static get darkGray() {
    return 0xFF838383;
  }
  static get blue() {
    return 0xFF0000FF;
  }
  static get green() {
    return 0xFF00FF00;
  }
  static get darkGreen() {
    return 0xFF006400;
  }
  static get black() {
    return 0xFF000000;
  }
  static get yellow() {
    return 0xFFFFFF00;
  }
  static get white() {
    return 0xFFFFFFFF;
  }
  static get red() {
    return 0xFFFF0000;
  }
  static get darkRed() {
    return 0xFF8b0000;
  }
  static get maroon() {
    return 0xFF800000;
  }
  static get olive() {
    return 0xFF808000;
  }
  static get orange() {
    return 0xFFFFA500;
  }
  static get darkOrange() {
    return 0xFFFF8C00;
  }
  static get lightBlue() {
    return 0xFFADD8E6;
  }
  static get darkBlue() {
    return 0xFF0000A0;
  }
  static get purple() {
    return 0xFF800080;
  }
  static get whitesmoke() {
    return 0xFFF5F5F5;
  }
  static get navy() {
    return 0xFF000080;
  }
  static get cyan() {
    return 0xFF00ffff;
  }
  static get filteredRows() {
    return 0xff1f77b4;
  }
  static get filteredOutRows() {
    return Color.lightLightGray;
  }
  static get selectedRows() {
    return Color.darkOrange;
  }
  static get missingValueRows() {
    return Color.filteredOutRows;
  }
  static get mouseOverRows() {
    return 0xFFAAAAAA;
  }
  static get currentRow() {
    return 0xFF38B738;
  }
  static get histogramBar() {
    return Color.filteredRows;
  }
  static get barChart() {
    return 0xFF24A221;
  }
  static get scatterPlotMarker() {
    return 0xFF40699c;
  }
  static get scatterPlotSelection() {
    return 0x80323232;
  }
  static get scatterPlotZoom() {
    return 0x80626200;
  }
  static get areaSelection() {
    return Color.lightBlue;
  }
  static get rowSelection() {
    return 0x60dcdca0;
  }
  static get colSelection() {
    return 0x60dcdca0;
  }
  static get areaZoom() {
    return 0x80323232;
  }
  static get gridWarningBackground() {
    return 0xFFFFB9A7;
  }
  static get success() {
    return 0xFF3cb173;
  }
  static get failure() {
    return 0xFFeb6767;
  }
}
/** Tree view node.
 * Sample: {@link https://public.datagrok.ai/js/samples/ui/tree-view}
 * */
export class TreeViewNode {
  /** @constructs {TreeView} */
  constructor(d) {
    this.d = d;
  }
  /** Creates new nodes tree
   * @returns {TreeViewNode} */
  static tree() {
    return toJs(api.grok_TreeViewNode_Tree());
  }
  /** Visual root.
   * @type {HTMLElement} */
  get root() {
    return api.grok_TreeViewNode_Root(this.d);
  }
  /** Caption label.
   * @type {HTMLElement} */
  get captionLabel() {
    return api.grok_TreeViewNode_CaptionLabel(this.d);
  }
  /** Check box.
   * @type {null|HTMLElement} */
  get checkBox() {
    return api.grok_TreeViewNode_CheckBox(this.d);
  }
  /** Returns 'true' if checked
   * @returns {boolean} */
  get checked() {
    return api.grok_TreeViewNode_Get_Checked(this.d);
  }
  set checked(checked) {
    api.grok_TreeViewNode_Set_Checked(this.d, checked);
  }
  /** Returns node text.
   * @returns {string} */
  get text() {
    return api.grok_TreeViewNode_Text(this.d);
  }
  /** Node value.
   * @type {Object}
   */
  get value() {
    return api.grok_TreeViewNode_Get_Value(this.d);
  }
  ;
  set value(v) {
    api.grok_TreeViewNode_Set_Value(this.d, v);
  }
  ;
  /** Gets all node items.
   * @returns {Array<TreeViewNode>} */
  get items() {
    return api.grok_TreeViewNode_Items(this.d).map((i) => toJs(i));
  }
  /** Add new group to node.
   * @param {string} text
   * @param {object} value
   * @param {boolean} expanded
   * @returns {TreeViewNode} */
  group(text, value = null, expanded = true) {
    return toJs(api.grok_TreeViewNode_Group(this.d, text, value, expanded));
  }
  /** Add new item to node.
   * @param {string} text
   * @param {object} value
   * @returns {TreeViewNode} */
  item(text, value = null) {
    return toJs(api.grok_TreeViewNode_Item(this.d, text, value));
  }
  /** Enables checkbox on node
   * @param {boolean} checked */
  enableCheckBox(checked = false) {
    api.grok_TreeViewNode_EnableCheckBox(this.d, checked);
  }
}
export class RangeSlider extends DartWidget {
  static create() {
    return toJs(api.grok_RangeSlider());
  }
  /** Gets minimum range value.
   * @type {number}
   */
  get minRange() {
    return api.grok_RangeSlider_Get_MinRange(this.d);
  }
  ;
  /** Gets maximum range value.
   * @type {number}
   */
  get maxRange() {
    return api.grok_RangeSlider_Get_MaxRange(this.d);
  }
  ;
  /** Gets minimum value.
   * @type {number}
   */
  get min() {
    return api.grok_RangeSlider_Get_Min(this.d);
  }
  ;
  /** Gets maximum value.
   * @type {number}
   */
  get max() {
    return api.grok_RangeSlider_Get_Max(this.d);
  }
  ;
  /** Sets values to range slider.
   * @param {number} minRange
   * @param {number} maxRange
   * @param {number} min
   * @param {number} max */
  setValues(minRange, maxRange, min, max) {
    api.grok_RangeSlider_SetValues(this.d, minRange, maxRange, min, max);
  }
  ;
  /** @returns {Observable} */
  get onValuesChanged() {
    return observeStream(api.grok_RangeSlider_Get_OnValuesChanged(this.d));
  }
}
export class HtmlTable extends DartWidget {
  /** @constructs {HtmlTable} */
  constructor(d) {
    super(d);
  }
  /** Creates a visual table based on [items], [renderer], and [columnNames]
   * @returns {HtmlTable} */
  static create(items, renderer, columnNames = null) {
    return toJs(api.grok_HtmlTable(items, renderer !== null ? (object, ind) => renderer(toJs(object), ind) : null, columnNames));
  }
  /** Removes item */
  remove(item) {
    api.grok_HtmlTable_Remove(this.d, item);
  }
}
export class ColumnComboBox extends DartWidget {
  /** @constructs {ColumnComboBox} */
  constructor(d) {
    super(d);
  }
  /** Creates a column combo box with specified [dataframe] and [predicate]
   * @returns {ColumnComboBox} */
  static create(dataframe, predicate) {
    return toJs(api.grok_ColumnComboBox(dataframe.d, (x) => predicate(toJs(x))));
  }
  /** @type {boolean} */
  get vertical() {
    return api.grok_ColumnComboBox_Get_Vertical(this.d);
  }
  /** @type {string} */
  get caption() {
    return api.grok_ColumnComboBox_Get_Caption(this.d);
  }
  set caption(c) {
    api.grok_ColumnComboBox_Set_Caption(this.d, c);
  }
  /** @type {Property} */
  get property() {
    return toJs(api.grok_ColumnComboBox_Get_Property(this.d));
  }
}
