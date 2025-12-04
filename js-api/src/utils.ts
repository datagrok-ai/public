import {Balloon} from './widgets';
import {Color} from './color';
import {toJs} from './wrappers';
import {COLUMN_TYPE, MARKER_TYPE, ViewerType, ColumnType} from './const';
import type {Point, Rect} from './grid';
import {IDartApi} from './api/grok_api.g';
import * as rxjs from 'rxjs';
import {StreamSubscription} from './events';
import {Column, DataFrame} from './dataframe';
import {TableView} from './views/view';
import wu from 'wu';
import { MapProxy } from './proxies';
import {ColumnInfo, DataConnection, TableInfo} from "./entities";

declare let DG: any;
declare let grok: any;

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

declare global {
  interface CanvasRenderingContext2D {
    setFillStyle(fill: string | CanvasGradient | CanvasPattern): CanvasRenderingContext2D;

    setStrokeStyle(stroke: string | CanvasGradient | CanvasPattern): CanvasRenderingContext2D;

    roundRect(x: number, y: number, w: number, h: number, r: number): CanvasRenderingContext2D

    line(x1: number, y1: number, x2: number, y2: number, color: number): CanvasRenderingContext2D;

    lines(points: Iterable<Point>): CanvasRenderingContext2D;

    /**
     * Use stroke() or fill() after.
     */
    polygon(pa: Point[]): CanvasRenderingContext2D;
  }

  interface HTMLCanvasElement {
    g2(): CanvasRenderingContext2D;
  }
}

if (typeof HTMLCanvasElement != 'undefined') {
  HTMLCanvasElement.prototype.g2 = function () {
    return this.getContext('2d')!;
  }

  CanvasRenderingContext2D.prototype.setFillStyle = function (fill: string | CanvasGradient | CanvasPattern) {
    this.fillStyle = fill;
    return this;
  }

  CanvasRenderingContext2D.prototype.setStrokeStyle = function (stroke: string | CanvasGradient | CanvasPattern) {
    this.strokeStyle = stroke;
    return this;
  }

  CanvasRenderingContext2D.prototype.roundRect = function (x: number, y: number, w: number, h: number, r: number) {
    if (w < 2 * r) r = w / 2;
    if (h < 2 * r) r = h / 2;
    this.beginPath();
    this.moveTo(x + r, y);
    this.arcTo(x + w, y, x + w, y + h, r);
    this.arcTo(x + w, y + h, x, y + h, r);
    this.arcTo(x, y + h, x, y, r);
    this.arcTo(x, y, x + w, y, r);
    this.closePath();
    return this;
  }

  CanvasRenderingContext2D.prototype.line = function (x1, y1, x2, y2, color) {
    this.beginPath();
    this.strokeStyle = Color.toRgb(color);
    this.moveTo(x1, y1);
    this.lineTo(x2, y2);
    this.stroke();

    return this;
  }

  CanvasRenderingContext2D.prototype.lines = function (points: Iterable<Point>) {
    let first = true;
    for (const point of points) {
      if (first) {
        this.moveTo(point.x, point.y);
        first = false;
      } else {
        this.lineTo(point.x, point.y);
      }
    }
    return this;
  }

  CanvasRenderingContext2D.prototype.polygon = function (pa: Point[]) {
    this.beginPath();

    const last_p = pa[pa.length - 1]
    this.moveTo(last_p.x, last_p.y);

    for (let p of pa) {
      this.lineTo(p.x, p.y);
    }
    this.closePath();

    return this;
  }
}


export namespace Paint {

  export function coordinateGrid(g: CanvasRenderingContext2D, worldBounds: Rect, xAxis?: Rect, yAxis?: Rect, viewport?: Rect) {
    if (xAxis != null)
      horzAxis(g, worldBounds.minX, worldBounds.maxX, xAxis.x, xAxis.y, xAxis.width, xAxis.height);
    if (yAxis != null)
      vertAxis(g, worldBounds.minY, worldBounds.maxY, yAxis.x, yAxis.y, yAxis.width, yAxis.height);
  }

  export function horzAxis(g: CanvasRenderingContext2D, min: number, max: number, x: number, y: number, w: number, h: number, log: boolean = false, inverse: boolean = false) {
    api.grok_Paint_HorzAxis(g, min, max, x, y, w, h, log, inverse);
  }

  export function vertAxis(g: CanvasRenderingContext2D, min: number, max: number, x: number, y: number, w: number, h: number, log: boolean = false, inverse: boolean = false) {
    api.grok_Paint_VertAxis(g, min, max, x, y, w, h, log, inverse);
  }

  export function roundRect(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, r: number): CanvasRenderingContext2D {
    if (w < 2 * r) r = w / 2;
    if (h < 2 * r) r = h / 2;
    g.beginPath();
    g.moveTo(x + r, y);
    g.arcTo(x + w, y, x + w, y + h, r);
    g.arcTo(x + w, y + h, x, y + h, r);
    g.arcTo(x, y + h, x, y, r);
    g.arcTo(x, y, x + w, y, r);
    g.closePath();
    return g;
  }

  /** Renders a marker */
  export function marker(g: CanvasRenderingContext2D, markerType: MARKER_TYPE, x: number, y: number, color: number | string, size: number) {
    const c: number = typeof color === 'string' ? Color.fromHtml(color) : color;
    api.grok_Paint_Marker(g, markerType, x, y, c, size);
  }

  export function markerTypes(): string[] {
    return api.grok_Marker_Types();
  }

  export function markerTypeIndexes(): {[key: string]: number} {
    return new MapProxy(api.grok_Marker_Type_Indexes()) as unknown as {[key: string]: number};
  }

  /** Renders a PNG image from bytes */
  export function pngImage(g: CanvasRenderingContext2D, bounds: Rect, imageBytes: Uint8Array) {
    let r = new FileReader();
    r.readAsBinaryString(new Blob([new Uint8Array(imageBytes)]));

    r.onload = function() {
      let img = new Image();

      img.onload = function() {
        bounds = new DG.Rect(bounds.x, bounds.y, bounds.width, bounds.height).fit(img.width, img.height)
        g.drawImage(img, bounds.x, bounds.y, bounds.width, bounds.height);
      }

      img.src = 'data:image/jpeg;base64,' + window.btoa(r.result as string);
    };
  }
}

export class Utils {
  /** @param {Iterable} iterable*/
  static firstOrNull<T>(iterable: Iterable<T>): T | null {
    let first = iterable[Symbol.iterator]().next();
    return first.done ? null : first.value;
  }

  /** Returns an 'identity' array where the element in idx-th position is equals to idx. */
  static identity(length: number) {
    let res = new Uint32Array(length);
    for (let i = 0; i < length; i++)
      res[i] = i;
    return res;
  }

  /** Returns random element from the array. Useful for demo data, tests, etc. */
  static random<T>(items: T[]): T {
    return items[Math.floor(Math.random() * items.length)];
  }

  static replaceAll(string: string, search: string, replace: string) {
    return string.split(search).join(replace);
  }

  static isEmpty(s?: string) {
    return !s || s === '';
  }

  static nullIfEmpty(s?: string) {
    return Utils.isEmpty(s) ? null : s;
  }

  /** Shows "Open File" dialog, and lets you process the result. */
  static openFile(options: {accept: string, open: (file: File) => void}): void {
    let input: HTMLInputElement = document.createElement('input');
    input.type = 'file';
    input.multiple = false;
    input.accept = options.accept ?? '';
    input.onchange = _this => {
        options.open(input.files![0]);
    };
    input.click();
  }

  /** Shows "Open File" dialog, and lets you process the result. */
  static openFileBytes(options: {accept: string, open: (bytes: Uint8Array) => void}): void {
    Utils.openFile({
      accept: options.accept,
      open: (f) => {
        var reader = new FileReader();
        reader.onload = () => {
          options.open(new Uint8Array(reader.result as ArrayBuffer));
        };
        reader.readAsArrayBuffer(f);
      }
    })
  }

  /** Downloads the specified content locally */
  static download(filename: string, content: BlobPart, contentType?: string) {
    contentType = contentType ?? 'application/octet-stream';
    const a = document.createElement('a');
    const blob = new Blob([content], {'type':contentType});
    a.href = window.URL.createObjectURL(blob);
    a.download = filename;
    a.click();
  }

  /** Loads the specified common libraries, if they were not loaded already.
   * Use it in plugins when a big JS dependency is used infrequently,
   * such as exporting to Excel, or showing a 3d structure.
   * Example: `loadJsCss(['common/exceljs.min.js', 'common/exceljs.min.css'])` */
  static async loadJsCss(files: string[]): Promise<null> {
    return toJs(api.grok_Utils_LoadJsCss(files));
  }

  static streamToObservable<T = any>(dartStream: any): rxjs.Observable<T> {
    return rxjs.fromEventPattern(
      function (handler) {
        return api.grok_Stream_Listen(dartStream, function (x: any) {
          handler(x);
        });
      },
      function (handler, dart) {
        new StreamSubscription(dart).cancel();
      }
    );
  }

  static getJsonValueType(x: any): ColumnType | null {
    if (!x)
      return null;

    if (typeof x === 'string')
      return DG.TYPE.STRING;
    if (typeof x === 'number')
      return DG.TYPE.FLOAT;
    if (typeof x === 'boolean')
      return DG.TYPE.BOOL;
    throw 'Not supported type';
  }

  /** Expands a column of JSON strings into multiple columns, one per unique key */
  static jsonToColumns(column: Column<string>) {
    const rows = column.dataFrame.rowCount;
    const jsons = range(rows).map(i => column.isNone(i) ? null : JSON.parse(column.get(i)!)).toArray();

    for (let i = 0; i < rows; i++) {
      for (const [key, value] of Object.entries(jsons[i])) {
        if (value !== null && !column.dataFrame.columns.contains(key))
          column.dataFrame.columns
            //@ts-ignore
            .addNew(key, Utils.getJsonValueType(value)!)
            .init(j => jsons[j] === null ? null : jsons[j][key]);
      }
    }
  }

  static async executeTests(testsParams: { package: any, params: any }[], stopOnFail?: boolean): Promise<any> {
    console.log(`********** Entered executeTests func`);
    let failed = false;
    let csv = "";
    let verbosePassed = "";
    let verboseSkipped = "";
    let verboseFailed = "";
    let countPassed = 0;
    let countSkipped = 0;
    let countFailed = 0;
    let resultDF: DataFrame | undefined = undefined;
    let lastTest: any = null;
    let res: string  = '';
    try {
      for (let testParam of testsParams) {
        lastTest = testParam;
        let df: DataFrame = await grok.functions.call(testParam.package + ':test', testParam.params);
        let flakingCol = DG.Column.fromType(DG.COLUMN_TYPE.BOOL, 'flaking', df.rowCount);
        df.columns.add(flakingCol);
        let packageNameCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'package', Array(df.rowCount).fill(testParam.package));
        df.columns.add(packageNameCol);
        if (df.rowCount === 0) {
          verboseFailed += `Test result : Invocation Fail : ${testParam.params.category}: ${testParam.params.test}\n`;
          countFailed += 1;
          failed = true;
          continue;
        }

        let row = df.rows.get(0);
        if (df.rowCount > 1) {
          let unhandledErrorRow = df.rows.get(1);
          if (!unhandledErrorRow.get("success")) {
            unhandledErrorRow["category"] = row.get("category");
            unhandledErrorRow["name"] = row.get("name");
            row = unhandledErrorRow;
          }
        }
        const category = row.get("category");
        const testName = row.get("name");
        const time = row.get("ms");
        const result = row.get("result");
        const success = row.get("success");
        const skipped = row.get("skipped");
        row["flaking"] = success && DG.Test.isReproducing;

        df.changeColumnType('result', COLUMN_TYPE.STRING);
        df.changeColumnType('logs', COLUMN_TYPE.STRING);
        df.changeColumnType('memoryDelta', COLUMN_TYPE.BIG_INT);

        if (resultDF === undefined)
          resultDF = df;
        else
          resultDF = resultDF.append(df);

        if (row["skipped"]) {
          verboseSkipped += `Test result : Skipped : ${time} : ${category}: ${testName} :  ${result}\n`;
          countSkipped += 1;
        }
        else if (row["success"]) {
          verbosePassed += `Test result : Success : ${time} : ${category}: ${testName} :  ${result}\n`;
          countPassed += 1;
        }
        else {
          verboseFailed += `Test result : Failed : ${time} : ${category}: ${testName} :  ${result}\n`;
          countFailed += 1;
          failed = true;
        }
        if ((success !== true && skipped !== true) && stopOnFail)
          break;
      }
      if (DG.Test.isInDebug) {
        console.log('on browser closing debug point');
        debugger
      }
      res = Utils.createResultsCsv(resultDF);

    } catch (e) {
      failed = true;
      verboseFailed = lastTest ? `category: ${lastTest.params.category}, name: ${lastTest.params.test}, error: ${e}, ${await DG.Logger.translateStackTrace((e as any).stack)}` :
        `test: null, error: ${e}, ${await DG.Logger.translateStackTrace((e as any).stack)}`;
    }

    return {
      failed: failed,
      verbosePassed: verbosePassed,
      verboseSkipped: verboseSkipped,
      verboseFailed: verboseFailed,
      passedAmount: countPassed,
      skippedAmount: countSkipped,
      failedAmount: countFailed,
      csv: res,
      // df: resultDF?.toJson()
    };
  }
  
  /**
   * Detects natural hierarchy in the list of categorical columns based on their values.
   * Mind that the function can be n^2 in the worst case scenario, although in most cases, it is quite fast
   * as it exits early if missmatch is found.
   * @param columns - List of columns to detect hierarchy in.
   * @param maxDepth - Maximum depth of the hierarchy to detect.
   * @returns 
   */
  static detectColumnHierarchy(columns: Column[], maxDepth: number = 3): string[] {
    return api.grok_Utils_DetectColumnHierarchy(columns.map((c) => c.dart), maxDepth);
  }

  static createResultsCsv(resultDF?: DataFrame): string {
    if (resultDF) {
      const bs = DG.BitSet.create(resultDF.rowCount)
      bs.setAll(true);
      for (let i = 0; i < resultDF.rowCount; i++) {
        if (resultDF.rows.get(i).get('category') === 'Unhandled exceptions') {
          bs.set(i, false);
        }
      }
      resultDF = resultDF.clone(bs);
      return resultDF.toCsv();
    }
    return '';
  }

  static randomString(length: number): string {
    const chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
    let result = '';
    const charsLen = chars.length;

    for (let i = 0; i < length; i++)
      result += chars.charAt(Math.floor(Math.random() * charsLen));

    return result;
  }
}


/** Html-related utilities */
export namespace HtmlUtils {

  /** Renders the element to canvas */
  export function renderToCanvas(element: HTMLElement, ratio: number = 1): Promise<HTMLCanvasElement> {
    return api.grok_HtmlUtils_RenderToCanvas(element, ratio);
  }

  export function htmlGetBounds(element: HTMLElement): Rect {
    return DG.Rect.fromDart(api.grok_HtmlUtils_HtmlGetBounds(element));
  }
}

export function _isDartium() {
  if (typeof document.head == 'undefined')
    return false;
  return document.head.querySelectorAll('script[src*=".dart.js"]').length == 0;
}


function* _range(length: number): IterableIterator<number> {
  for (let i = 0; i < length; i++)
    yield i;
}

/** Generates [count] increasing integer numbers, starting with 0. */
export function range(length: number): wu.WuIterable<number> {

  return wu(_range(length));
}

/** Returns an 'identity' array where the element in idx-th position is equals to idx. */
export function identity(length: number) {
  let res = new Uint32Array(length);
  for (let i = 0; i < length; i++)
    res[i] = i;
  return res;
}

/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return api.grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/

/** Times the execution of function f
 * @param {string} name - a label for the execution time to display
 * @param {Function} f - function with no parameters that will get measured
 * @returns {value} - a value which f returns
 * */
export function time(name: string, f: Function) {
  let start = new Date();
  let result = f();
  let stop = new Date();
  // @ts-ignore
  console.log(`${name}: ${stop - start} ms`);
  // @ts-ignore
  new Balloon().info(`${name}: ${stop - start} ms`);
  return result;
}

/** Times the execution of asyncronous function f
 * @async
 * @param {string} name - a label for the execution time to display
 * @param {Function} f - async function with no parameters that will get measured
 * @returns {Promise<value>} - a promise for the value which f returns
 * */
export async function timeAsync(name: string, f: Function) {
  let start = new Date();
  let result = await f();
  let stop = new Date();
  // @ts-ignore
  console.log(`${name}: ${stop - start} ms`);
  // @ts-ignore
  new Balloon().info(`${name}: ${stop - start} ms`);
  return result;
}

export function _identityInt32(length: number): Int32Array {
  let values = new Int32Array(length);
  for (let i = 0; i < length; i++)
    values[i] = i;
  return values;
}


/**
 * Least recently used cache.
 * Inspired by https://github.com/Yomguithereal/mnemonist/blob/master/lru-cache.js
 * */
export class LruCache<K = any, V = any> {
  private capacity: number;
  public onItemEvicted: Function | null;
  private items: {};
  private tail: number;
  private forward: Uint16Array;
  private backward: Uint16Array;
  private V: any[];
  private K: any[];
  private size: number;
  private head: number;

  constructor(capacity: number = 100) {
    this.capacity = capacity;
    this.forward = new Uint16Array(this.capacity);
    this.backward = new Uint16Array(this.capacity);
    this.V = new Array(this.capacity);
    this.K = new Array(this.capacity);
    this.size = 0;
    this.head = 0;
    this.tail = 0;
    this.items = {};
    this.onItemEvicted = null;
  }

  /**
   * Splays a value on top.
   * @param {number} pointer - Pointer of the value to splay on top.
   * @return {LruCache}
   */
  splayOnTop(pointer: number): LruCache<K, V> {
    let oldHead = this.head;

    if (this.head === pointer)
      return this;

    let previous = this.backward[pointer];
    let next = this.forward[pointer];

    if (this.tail === pointer)
      this.tail = previous;
    else
      this.backward[next] = previous;

    this.forward[previous] = next;
    this.backward[oldHead] = pointer;
    this.head = pointer;
    this.forward[pointer] = oldHead;

    return this;
  }


  /**
   * Checks whether the key exists in the cache.
   *
   * @param  {any} key   - Key.
   */
  has(key: any): boolean {
    return key in this.items;
  }

  /**
   * Sets the value for the given key in the cache.
   *
   * @param  {any} key   - Key.
   * @param  {any} value - Value.
   */
  set(key: any, value: any): void {

    // The key already exists, we just need to update the value and splay on top
    // @ts-ignore
    let pointer = this.items[key];

    if (typeof pointer !== 'undefined') {
      this.splayOnTop(pointer);
      this.V[pointer] = value;

      return;
    }

    // The cache is not yet full
    if (this.size < this.capacity) {
      pointer = this.size++;
    }

    // Cache is full, we need to drop the last value
    else {
      pointer = this.tail;
      this.tail = this.backward[pointer];
      if (this.onItemEvicted != null)
        this.onItemEvicted(this.V[pointer]);
      // @ts-ignore
      delete this.items[this.K[pointer]];
    }

    // Storing key & value
    // @ts-ignore
    this.items[key] = pointer;
    this.K[pointer] = key;
    this.V[pointer] = value;

    // Moving the item at the front of the list
    this.forward[pointer] = this.head;
    this.backward[this.head] = pointer;
    this.head = pointer;
  }

  /**
   * Gets the value attached to the given key, and makes it the most recently used item.
   *
   * @param  {any} key - Key.
   */
  get(key: any): V | undefined {
    // @ts-ignore
    let pointer = this.items[key];

    if (typeof pointer === 'undefined')
      return undefined;

    this.splayOnTop(pointer);

    return this.V[pointer];
  }

  /**
   * Returns the value with the specified key, if it already exists in the cache,
   * or creates a new one by calling the provided function.
   *
   * @param  {any} key   - Key.
   * @param  {Function} createFromKey - Function to create a new item.
   * @return {any}
   */
  getOrCreate(key: K, createFromKey: (key: K) => V): V {
    let value = this.get(key);
    if (value !== undefined)
      return value;
    else {
      let item = createFromKey(key);
      this.set(key, item);
      return item;
    }
  }
}


export function _options(element: HTMLElement, options: any): HTMLElement {
  if (options == null)
    return element;
  if (typeof options === 'string')
    element.className += ` ${options.replace(/,/g, ' ')}`;
  if (options.id != null)
    element.id = options.id;
  if (options.classes != null)
    element.className += ` ${options.classes}`;
  if (options.style != null)
    Object.assign(element.style, options.style);
  if (options.processNode != null)
    options.processNode(element);
  if (options.onClick != null)
    element.addEventListener('click', options.onClick);
  return element;
}

export function format(x: number, format?: string): string {
  return api.grok_Utils_FormatNumber(x, format);
}

/* Waits [ms] milliseconds */
export async function delay(ms: number) {
  await new Promise((r) => setTimeout(r, ms));
}

/** Autotest-related helpers */
export namespace Test {

  /** Specifies whether the TestManager works in benchmark on unit test mode.
   * This is controlled by the TestManager's "Benchmark" check box.
   *
   * Use this variable in the individual tests to control whether it
   * runs as a quick unit test (for instance, set small number of iterations)
   * or a benchmark (set big number of iterations).
   *
   * Results of the benchmarks, including time, are captured during the overnight
   * automated testing and are reported to the metrics database. This allows to
   * identify performance regressions, compare how fast tests perform under
   * different conditions, etc.
   * */
  export let isInBenchmark = false;
  export let isReproducing = false;
  export let isInDebug = false;
  export let isCiCd = false;
  export let isProfiling = false;

  export function getTestDataGeneratorByType(type: string) {
    return api.grok_Test_GetTestDataGeneratorByType(type);
  }

  export function getInputTestDataGeneratorByType(inputType: string) {
    return api.grok_Test_GetInputTestDataGeneratorByType(inputType);
  }

  export async function testViewerProperties(df: DataFrame, viewerType: ViewerType, properties: { [key: string]: any}): Promise<TableView> {
    const tv = TableView.create(df, true);
    const viewer = tv.addViewer(viewerType);
    viewer.setOptions(properties);
    await delay(1000);
    return tv;
  }

  export async function testViewer(viewerType: ViewerType) {
    await api.grok_Test_RunViewerTest(viewerType);
  }
}

// Completer class based on https://stackoverflow.com/a/67007151
export class Completer<T> {
	public readonly promise: Promise<T>;
	public complete: (value: (PromiseLike<T> | T)) => void;
	public reject: (reason?: any) => void;

	public constructor() {
    this.complete = (_) => {};
    this.reject = (_) => {};
    this.promise = new Promise<T>((resolve, reject) => {
        this.complete = resolve;
        this.reject = reject;
    })
	}
}

/**
 * Represents metadata for a specific database schema in Datagrok.
 *
 * `DbSchemaInfo` provides access to schema-level information such as tables,
 * user and LLM comments, and utilities for annotating tables and columns.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbSchemaInfo {
  public dart: any;

  /**
   * Creates a new `DbSchemaInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a database schema.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * Schema name (e.g., `"public"`, `"dbo"`).
   *
   * @returns The name of this schema.
   */
  get name(): string {
    return api.grok_DbSchemaInfo_Get_Name(this.dart);
  }

  /**
   * @returns A comment string, or an empty string if none exists.
   */
  get comment(): string {
    return api.grok_DbSchemaInfo_Get_Comment(this.dart);
  }

  /**
   * @returns User-defined comment used by AI query builder.
   */
  get llmComment(): string {
    return api.grok_DbSchemaInfo_Get_LlmComment(this.dart);
  }

  /**
   * The parent `DataConnection` this schema belongs to.
   *
   * @returns A {@link DataConnection} instance.
   */
  get connection(): DataConnection {
    return api.grok_DbSchemaInfo_Get_Connection(this.dart);
  }

  /**
   * Lists all tables belonging to this schema.
   *
   * @returns A promise resolving to an array of {@link TableInfo} objects.
   */
  get tables(): Promise<TableInfo[]> {
    return api.grok_DbSchemaInfo_Get_Tables(this.dart);
  }

  // ─────────────────────────────── SETTERS ───────────────────────────────

  setComment(comment: string): Promise<void> {
    return api.grok_DbSchemaInfo_Set_Comment(this.dart, comment);
  }

  setLlmComment(llmComment: string): Promise<void> {
    return api.grok_DbSchemaInfo_Set_LlmComment(this.dart, llmComment);
  }

  // ─────────────────────────────── ANNOTATION METHODS ───────────────────────────────

  /**
   * Adds or updates metadata for a specific table in this schema.
   * It's just only meta annotation which is stored in Datagrok and real database is not affected.
   *
   * @param table - Either a {@link TableInfo} instance or a table name.
   * @param props - Table-level properties such as comments, row counts, and domains.
   *
   * @returns A promise that resolves when the annotation is saved.
   *
   * @example
   * ```ts
   * await schema.annotateTable("customers", {
   *   comment: "Master table containing customer profiles",
   *   rowCount: 120345,
   * });
   * ```
   */
  annotateTable(
      table: TableInfo | string,
      props: DbTableProperties,
  ): Promise<void> {
    return api.grok_DbSchemaInfo_AnnotateTable(
        this.dart,
        table instanceof TableInfo ? table.friendlyName : table,
        props,
    );
  }

  /**
   * Adds or updates metadata for a column within a table.
   * It's just only meta annotation which is stored in Datagrok and real database is not affected.
   *
   * @param table - A {@link TableInfo} instance or table name.
   * @param column - A {@link ColumnInfo} instance or column name.
   * @param props - Column-level annotations such as uniqueness, min/max values, and comments.
   *
   * @returns A promise that resolves when the annotation is saved.
   *
   * @example
   * ```ts
   * await schema.annotateColumn("orders", "order_id", {
   *   isUnique: true,
   *   comment: "Primary key for the orders table",
   * });
   * ```
   */
  annotateColumn(
      table: TableInfo | string,
      column: ColumnInfo | string,
      props: DbColumnProperties,
  ): Promise<void> {
    return api.grok_DbSchemaInfo_AnnotateColumn(
        this.dart,
        table instanceof TableInfo ? table.friendlyName : table, //  table.friendlyName corresponds to real db table name
        column instanceof ColumnInfo ? column.name : column, // column.name corresponds to real db col name
        props,
    );
  }
}


/**
 * Represents metadata for a database relation in Datagrok.
 *
 * `DbRelationInfo` describes how two tables are connected, including the
 * schema/table pairs, participating columns, cardinality, and whether this
 * relation is considered a primary navigation path.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbRelationInfo implements DbRelationProperties {
  public dart: any;

  /**
   * Creates a new `DbRelationInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a relation.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * The name of this relation.
   */
  get name(): string {
    return api.grok_DbRelationInfo_Get_Name(this.dart);
  }

  /**
   * User-defined comment describing this relation.
   */
  get comment(): string {
    return api.grok_DbRelationInfo_Get_Comment(this.dart);
  }

  /**
   * LLM-generated annotation or AI summary for this relation.
   */
  get llmComment(): string {
    return api.grok_DbRelationInfo_Get_LlmComment(this.dart);
  }

  /**
   * The relation cardinality (e.g., `"one-to-many"`, `"many-to-one"`).
   */
  get cardinality(): string {
    return api.grok_DbRelationInfo_Get_Cardinality(this.dart);
  }

  /**
   * Indicates whether this relation is treated as a primary navigation path
   * (e.g., preferred for automated join suggestions).
   */
  get isPrimaryPath(): boolean {
    return toJs(api.grok_DbRelationInfo_Get_IsPrimaryPath(this.dart));
  }

  /**
   * Connection object this relation belongs to.
   */
  get connection(): DataConnection {
    return toJs(api.grok_DbRelationInfo_Get_Connection(this.dart));
  }

  /**
   * Schema name of the source table of this relation.
   */
  get fromSchema(): string {
    return api.grok_DbRelationInfo_Get_FromSchema(this.dart);
  }

  /**
   * Source table name of this relation.
   */
  get fromTable(): string {
    return api.grok_DbRelationInfo_Get_FromTable(this.dart);
  }

  /**
   * List of source table columns participating in the relation.
   */
  get fromColumns(): string[] {
    return api.grok_DbRelationInfo_Get_FromColumns(this.dart);
  }

  /**
   * Schema name of the target table of this relation.
   */
  get toSchema(): string {
    return api.grok_DbRelationInfo_Get_ToSchema(this.dart);
  }

  /**
   * Target table name of this relation.
   */
  get toTable(): string {
    return api.grok_DbRelationInfo_Get_ToTable(this.dart);
  }

  /**
   * List of target table columns participating in the relation.
   */
  get toColumns(): string[] {
    return api.grok_DbRelationInfo_Get_ToColumns(this.dart);
  }

  // ─────────────────────────────── MODIFIER METHODS ───────────────────────────────

  setComment(comment: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_Comment(this.dart, comment);
  }

  setLlmComment(llmComment: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_LlmComment(this.dart, llmComment);
  }

  setCardinality(cardinality: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_Cardinality(this.dart, cardinality);
  }

  setIsPrimaryPath(isPrimaryPath: boolean): Promise<void> {
    return api.grok_DbRelationInfo_Set_IsPrimaryPath(this.dart, isPrimaryPath);
  }

  setFromSchema(schema: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_FromSchema(this.dart, schema);
  }

  setFromTable(table: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_FromTable(this.dart, table);
  }

  setFromColumns(cols: string[]): Promise<void> {
    return api.grok_DbRelationInfo_Set_FromColumns(this.dart, cols);
  }

  setToSchema(schema: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_ToSchema(this.dart, schema);
  }

  setToTable(table: string): Promise<void> {
    return api.grok_DbRelationInfo_Set_ToTable(this.dart, table);
  }

  setToColumns(cols: string[]): Promise<void> {
    return api.grok_DbRelationInfo_Set_ToColumns(this.dart, cols);
  }
}


export interface DbRelationProperties {
  comment?: string;
  llmComment?: string;
  cardinality?: string;
  isPrimaryPath?: boolean;
  fromSchema?: string;
  fromTable?: string;
  fromColumns?: string[];
  toSchema?: string;
  toTable?: string;
  toColumns?: string[];
}

export interface DbColumnProperties {
  comment?: string;
  llmComment?: string;
  isUnique?: boolean;
  min?: number;
  max?: number;
  values?: string;
  sampleValues?: string;
  uniqueCount?: number;
}

export interface DbTableProperties {
  comment?: string;
  llmComment?: string;
  domains?: string;
  rowCount?: number;
}

/**
 * Represents metadata for a database connection in Datagrok.
 *
 * Provides access to high-level database metadata such as schemas,
 * relations, comments, and LLM-generated annotations. It also allows updating
 * connection-level annotations and adding new relation metadata.
 *
 * Mutation methods do not modify underlying database, they just create metadata in Datagrok.
 */
export class DbInfo {
  public dart: any;

  /**
   * Creates a new `DbInfo` wrapper.
   *
   * @param dart - The underlying Dart object representing a database connection.
   */
  constructor(dart: any) {
    this.dart = dart;
  }

  // ─────────────────────────────── GETTERS ───────────────────────────────

  /**
   * The name of the database. If database name is not set in connection details then this field is undefined.
   *
   * @returns The database name.
   */
  get name(): string | undefined {
    return api.grok_DbInfo_Get_Name(this.dart);
  }

  /**
   * @returns Connection comment text, or an empty string if none exists.
   */
  get comment(): string {
    return api.grok_DbInfo_Get_Comment(this.dart);
  }

  /**
   * @returns User-defined comment used by AI query builder.
   */
  get llmComment(): string {
    return api.grok_DbInfo_Get_LlmComment(this.dart);
  }

  /**
   * Lists all schemas available for this database connection.
   *
   * @returns A promise resolving to an array of {@link DbSchemaInfo} objects.
   */
  get schemas(): Promise<DbSchemaInfo[]> {
    return api.grok_DbInfo_Get_Schemas(this.dart);
  }

  /**
   * Lists all known relations for this connection.
   *
   * @returns A promise resolving to an array of {@link DbRelationInfo} objects.
   */
  get relations(): Promise<DbRelationInfo[]> {
    return api.grok_DbInfo_Get_Relations(this.dart);
  }

  /**
   * The underlying `DataConnection` object.
   *
   * @returns A {@link DataConnection} wrapper.
   */
  get connection(): DataConnection {
    return toJs(api.grok_DbInfo_Get_Connection(this.dart));
  }

  // ─────────────────────────────── SETTERS ───────────────────────────────

  setComment(comment: string): Promise<void> {
    return api.grok_DbInfo_Set_Comment(this.dart, comment);
  }

  setLlmComment(comment: string): Promise<void> {
    return api.grok_DbInfo_Set_LlmComment(this.dart, comment);
  }

  // ─────────────────────────────── MUTATION ───────────────────────────────

  /**
   * Creates and registers a new database relation for this connection.
   * It's just only meta annotation which is stored in Datagrok and real database is not affected.
   *
   * @param name - Relation name.
   * @param props - Relation metadata, including table/column mappings.
   * @returns A promise resolving to the newly created {@link DbRelationInfo}.
   *
   * @example
   * ```ts
   * const rel = await dbInfo.addRelation('fk_orders_customers', {
   *   cardinality: 'many-to-one',
   *   fromTable: 'orders',
   *   fromColumns: ['customer_id'],
   *   toTable: 'customers',
   *   toColumns: ['id'],
   * });
   * ```
   */
  addRelation(name: string, props: DbRelationProperties): Promise<DbRelationInfo> {
    return api.grok_DbInfo_AddRelation(this.dart, name, props);
  }

  /**
   * Removes all the meta data associated with the connection in Datagrok.
   */
  clearProperties(): Promise<void> {
    return api.grok_DbInfo_ClearProperties(this.dart);
  }
}
