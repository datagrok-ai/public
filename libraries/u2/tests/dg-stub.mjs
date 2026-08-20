/* Module hooks that serve a minimal `datagrok-api` to node, so the dg modules that need the
   platform at runtime can be tested headlessly. The real package does not load in node (its dayjs
   plugin imports carry no extension), and the platform base classes call into dart in their
   constructors.

   The stub is the platform's own contract, nothing more: the `DG.TYPE` strings copied from
   js-api's `const.ts`, the `JsInputBase` surface `DartInput` uses (root, caption, addValidator,
   fireInput/fireChanged, and `property`, which throws off a dart handle when nothing is bound),
   the `FileInfo` the file inputs build and read, and the `grok.dapi.files` / `grok.shell` calls
   they make. Everything on `dapi`, `shell` and `ui` is a mutable field a test replaces with its
   own function — import the stub module in the test and assign.

   Register it with `register('./dg-stub.mjs', import.meta.url)` before importing the module under
   test. */

const STUB = `
export const TYPE = {
  STRING: 'string', INT: 'int', FLOAT: 'double', NUM: 'num', BOOL: 'bool', DATE_TIME: 'datetime',
  BIG_INT: 'bigint', QNUM: 'qnum', OBJECT: 'object', FILE: 'file',
};

/** Enough of the entity for the file inputs: what they build (fromBytes) and what they read.
 * A file in a share carries the connection it lives on; one built out of local bytes does not. */
export class FileInfo {
  constructor(path, data, connection = null) {
    this.fullPath = path;
    this.path = path;
    this.fileName = path.substring(path.lastIndexOf('/') + 1);
    this.data = data;
    this.connection = connection;
    this.isFile = true;
    this.isDirectory = false;
  }

  static fromBytes(path, data) { return new FileInfo(path, data); }
}

/** The platform property the widget host mints out of a u2 PropertyLike: what it is created with,
 * plus the options record it is refined with. */
export class Property {
  static create(name, type, get, set, defaultValue) {
    const p = new Property();
    p.name = name;
    p.type = type;
    p.propertyType = type;
    p.get = get;
    p.set = set;
    p.defaultValue = defaultValue;
    return p;
  }

  fromOptions(options) {
    Object.assign(this, options);
    return this;
  }
}

/** The base widget contract the host subclasses: a root, the kill channel, and the introspection
 * members it overrides (js-api widgets/base.ts). */
export class Widget {
  constructor(root) {
    this._root = root;
    this.subs = [];
    this.temp = {};
    this._properties = [];
    this._functions = [];
    this._aiDescription = null;
    this.isDetached = false;
  }

  static fromRoot(root) { return new Widget(root); }

  get type() { return 'Unknown'; }

  get root() { return this._root; }

  getProperties() { return this._properties; }

  getFunctions() { return this._functions; }

  get aiDescription() { return this._aiDescription; }

  set aiDescription(x) { this._aiDescription = x; }

  onEvent(_eventId = null) { return {subscribe: () => ({unsubscribe() {}})}; }

  getWidgetStatus() {
    return {parts: {}, hitAreas: {}, shortcuts: {}, events: [], description: null, error: null};
  }

  toDart() {
    this.dart ??= {widget: this};
    return this.dart;
  }

  detach() {
    for (const s of this.subs)
      s.unsubscribe();
    this.isDetached = true;
  }
}

/** The base a handler subclasses (js-api ui.ts:1687): what it must override throws, what it may
 * override has the platform's default, and the registry is the list \`register\` appends to. */
export class ObjectHandler {
  static registered = [];

  get type() { throw 'Not defined.'; }

  get name() { return \`\${this.type} handler\`; }

  isApplicable(_x) { throw 'Not defined.'; }

  getCaption(x) { return \`\${x}\`; }

  static register(handler) { ObjectHandler.registered.push(handler); }
}

export class JsInputBase {
  constructor() {
    this.root = document.createElement('div');
    this.root.classList.add('ui-input-root');
    // validators live dart-side, on the proxy \`grok_InputBase_FromJS\` builds — a second proxy
    // (\`JsInputProxy.fromFunc\`) re-homes the input and carries its own list
    this.dart = {validators: []};
    this.fired = {input: 0, changed: 0};
  }

  /** The platform reads it off the dart handle, which throws when the proxy bound nothing. */
  get property() {
    if (this._property === undefined)
      throw new Error('no dart handle');
    return this._property;
  }

  set property(x) { this._property = x; }

  addValidator(v) { this.dart.validators.push(v); }

  fireInput() { this.fired.input++; }

  fireChanged() { this.fired.changed++; }
}
`;

const GROK_STUB = `
export const dapi = {
  files: {
    exists: async () => false,
    list: async () => [],
    readAsText: async () => '',
  },
};

export const shell = {
  info: () => {},
  warning: () => {},
  error: () => {},
  /** The current object — what a selection writes to, and what the context panel renders. */
  o: null,
  windows: {showContextPanel: false},
};

/** Ad-hoc function registration (js-api functions.ts:137): every registration, newest last, and
 * the Func it hands back — the platform's own is a dart wrapper over the same record. */
export const functions = {
  registrations: [],
  register(data) {
    functions.registrations.push(data);
    return {...data, name: data.signature.split(/[ (]/)[1]};
  },
};
`;

const UI_STUB = `
/** Every registration, newest last — a test drives the platform drop channel through it. */
export const drops = [];

export function makeDroppable(element, options) {
  drops.push({element, ...options});
}
`;

const URLS = {
  'datagrok-api/dg': ['u2-test:datagrok-api-dg', STUB],
  'datagrok-api/grok': ['u2-test:datagrok-api-grok', GROK_STUB],
  'datagrok-api/ui': ['u2-test:datagrok-api-ui', UI_STUB],
};
const SOURCES = new Map(Object.values(URLS));

export async function resolve(specifier, context, next) {
  const stub = URLS[specifier];
  if (stub)
    return {url: stub[0], format: 'module', shortCircuit: true};
  return next(specifier, context);
}

export async function load(url, context, next) {
  const source = SOURCES.get(url);
  if (source)
    return {source, format: 'module', shortCircuit: true};
  return next(url, context);
}
