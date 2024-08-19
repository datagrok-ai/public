/// this file was generated automatically from d4 classes declarations
import { toDart } from "../wrappers";
let api = <any>window;

export class UsageType {
  static CLICK = 'click';

  static MENU_CLICK = 'menu click';

  static DIALOG_OK = 'dialog ok';

  static DIALOG_SHOW = 'dialog show';

  static USER_REPORT_POSTED = 'user report posted';

}
export class ViewerEvent {
  public dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }
  static create(): ViewerEvent {
    return new ViewerEvent(api.grok_ViewerEvent_Create());
  }
  get viewer(): any { return api.grok_ViewerEvent_Get_viewer(this.dart); };
  set viewer(x: any) {api.grok_ViewerEvent_Set_viewer(this.dart, toDart(x)); }
  get type(): string { return api.grok_ViewerEvent_Get_type(this.dart); };
  set type(x: string) {api.grok_ViewerEvent_Set_type(this.dart, toDart(x)); }
  get eventFlag(): boolean { return api.grok_ViewerEvent_Get_eventFlag(this.dart); };
  set eventFlag(x: boolean) {api.grok_ViewerEvent_Set_eventFlag(this.dart, toDart(x)); }
  get filters(): {[index: string]: any} { return api.grok_ViewerEvent_Get_filters(this.dart); };
  set filters(x: {[index: string]: any}) {api.grok_ViewerEvent_Set_filters(this.dart, toDart(x)); }
  get row(): number { return api.grok_ViewerEvent_Get_row(this.dart); };
  set row(x: number) {api.grok_ViewerEvent_Set_row(this.dart, toDart(x)); }
  get mouseEvent(): any { return api.grok_ViewerEvent_Get_mouseEvent(this.dart); };
  set mouseEvent(x: any) {api.grok_ViewerEvent_Set_mouseEvent(this.dart, toDart(x)); }
  get bitset(): any { return api.grok_ViewerEvent_Get_bitset(this.dart); };

}
export class InputType {
  static Int = 'Int';

  static BigInt = 'BigInt';

  static Float = 'Float';

  static QNum = 'QNum';

  static Slider = 'Slider';

  static Bool = 'Bool';

  static Switch = 'Switch';

  static Text = 'Text';

  static TextArea = 'TextArea';

  static Markdown = 'Markdown';

  static Code = 'Code';

  static Search = 'Search';

  static Date = 'Date';

  static Map = 'Map';

  static File = 'File';

  static List = 'List';

  static Color = 'Color';

  static Column = 'Column';

  static Columns = 'Columns';

  static ColumnsMap = 'ColumnsMap';

  static Radio = 'Radio';

  static Choice = 'Choice';

  static MultiChoice = 'MultiChoice';

  static Table = 'Table';

  static Molecule = 'Molecule';

  static User = 'User';

  static UserGroups = 'UserGroups';

  static Dynamic = 'Dynamic';

  static Image = 'Image';

  static JsInputProxy = 'JsInputProxy';

}
export class GridCellStyleEx {
  public dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }
  static create(): GridCellStyleEx {
    return new GridCellStyleEx(api.grok_GridCellStyleEx_Create());
  }
  static DATA_TYPE = 'data type';

  static CELL_TYPE = 'cell type';

  static VERT_ALIGN_TOP = 'top';

  static VERT_ALIGN_CENTER = 'center';

  static VERT_ALIGN_BOTTOM = 'bottom';

  static HORZ_ALIGN_LEFT = 'left';

  static HORZ_ALIGN_CENTER = 'center';

  static HORZ_ALIGN_RIGHT = 'right';

  static STYLE_DEFAULT = 'default';

  static STYLE_TEXT = 'text';

  static STYLE_NUMBER = 'number';

  static STYLE_ROW_HEADER = 'row header';

  static STYLE_COL_HEADER = 'column header';

  static get defaultStyle(): any { return api.grok_GridCellStyle_Get_defaultStyle(); };
  static set defaultStyle(x: any) {api.grok_GridCellStyle_Set_defaultStyle(toDart(x)); }
  static get textStyle(): any { return api.grok_GridCellStyle_Get_textStyle(); };
  static set textStyle(x: any) {api.grok_GridCellStyle_Set_textStyle(toDart(x)); }
  static get numberStyle(): any { return api.grok_GridCellStyle_Get_numberStyle(); };
  static set numberStyle(x: any) {api.grok_GridCellStyle_Set_numberStyle(toDart(x)); }
  static get styles(): {[index: string]: any} { return api.grok_GridCellStyle_Get_styles(); };
  static set styles(x: {[index: string]: any}) {api.grok_GridCellStyle_Set_styles(toDart(x)); }
  get font(): string { return api.grok_GridCellStyle_Get_font(this.dart); };
  set font(x: string) {api.grok_GridCellStyle_Set_font(this.dart, toDart(x)); }
  get horzAlign(): string { return api.grok_GridCellStyle_Get_horzAlign(this.dart); };
  set horzAlign(x: string) {api.grok_GridCellStyle_Set_horzAlign(this.dart, toDart(x)); }
  get vertAlign(): string { return api.grok_GridCellStyle_Get_vertAlign(this.dart); };
  set vertAlign(x: string) {api.grok_GridCellStyle_Set_vertAlign(this.dart, toDart(x)); }
  /// When defined, overrides the default cell tooltip
  get tooltip(): string { return api.grok_GridCellStyle_Get_tooltip(this.dart); };
  set tooltip(x: string) {api.grok_GridCellStyle_Set_tooltip(this.dart, toDart(x)); }
  get cursor(): string { return api.grok_GridCellStyle_Get_cursor(this.dart); };
  set cursor(x: string) {api.grok_GridCellStyle_Set_cursor(this.dart, toDart(x)); }
  get textWrap(): string { return api.grok_GridCellStyle_Get_textWrap(this.dart); };
  set textWrap(x: string) {api.grok_GridCellStyle_Set_textWrap(this.dart, toDart(x)); }
  /// Marker to be shown when the value does not fit in the cell
  get marker(): string { return api.grok_GridCellStyle_Get_marker(this.dart); };
  set marker(x: string) {api.grok_GridCellStyle_Set_marker(this.dart, toDart(x)); }
  get textColor(): number { return api.grok_GridCellStyle_Get_textColor(this.dart); };
  set textColor(x: number) {api.grok_GridCellStyle_Set_textColor(this.dart, toDart(x)); }
  get backColor(): number { return api.grok_GridCellStyle_Get_backColor(this.dart); };
  set backColor(x: number) {api.grok_GridCellStyle_Set_backColor(this.dart, toDart(x)); }
  get marginLeft(): number { return api.grok_GridCellStyle_Get_marginLeft(this.dart); };
  set marginLeft(x: number) {api.grok_GridCellStyle_Set_marginLeft(this.dart, toDart(x)); }
  get marginRight(): number { return api.grok_GridCellStyle_Get_marginRight(this.dart); };
  set marginRight(x: number) {api.grok_GridCellStyle_Set_marginRight(this.dart, toDart(x)); }
  get marginTop(): number { return api.grok_GridCellStyle_Get_marginTop(this.dart); };
  set marginTop(x: number) {api.grok_GridCellStyle_Set_marginTop(this.dart, toDart(x)); }
  get marginBottom(): number { return api.grok_GridCellStyle_Get_marginBottom(this.dart); };
  set marginBottom(x: number) {api.grok_GridCellStyle_Set_marginBottom(this.dart, toDart(x)); }
  get textVertical(): boolean { return api.grok_GridCellStyle_Get_textVertical(this.dart); };
  set textVertical(x: boolean) {api.grok_GridCellStyle_Set_textVertical(this.dart, toDart(x)); }
  /// Applies to image columns only
  get imageScale(): number { return api.grok_GridCellStyle_Get_imageScale(this.dart); };
  set imageScale(x: number) {api.grok_GridCellStyle_Set_imageScale(this.dart, toDart(x)); }
  /// Applies to image columns only
  get opacity(): number { return api.grok_GridCellStyle_Get_opacity(this.dart); };
  set opacity(x: number) {api.grok_GridCellStyle_Set_opacity(this.dart, toDart(x)); }
  get clip(): boolean { return api.grok_GridCellStyle_Get_clip(this.dart); };
  set clip(x: boolean) {api.grok_GridCellStyle_Set_clip(this.dart, toDart(x)); }
  /// For 'html' cell types only
  get element(): any { return api.grok_GridCellStyle_Get_element(this.dart); };
  set element(x: any) {api.grok_GridCellStyle_Set_element(this.dart, toDart(x)); }
  /// When defined, the cell editor becomes a combo box with the specified values
  get choices(): Array<string> { return api.grok_GridCellStyle_Get_choices(this.dart); };
  set choices(x: Array<string>) {api.grok_GridCellStyle_Set_choices(this.dart, toDart(x)); }
}
export function renderMultipleHistograms(g: CanvasRenderingContext2D, bounds: any, histograms: Array<Int32List>, options?: {categoryColumn?: any, colors?: Array<number>, tension?: number, normalize?: boolean, markerSize?: number, fill?: boolean, minBin?: number, maxBin?: number, localMaximum?: boolean, highlightedHistogram?: number}): any
  { return api.grok_renderMultipleHistograms(toDart(g), toDart(bounds), toDart(histograms), toDart(options?.categoryColumn), toDart(options?.colors), toDart(options?.tension), toDart(options?.normalize), toDart(options?.markerSize), toDart(options?.fill), toDart(options?.minBin), toDart(options?.maxBin), toDart(options?.localMaximum), toDart(options?.highlightedHistogram)); }

