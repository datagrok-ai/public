/**
 * Chem's Crux module, recovering: a user is never asked to reload, a failure is handled inside, and what fails is
 * said. A wasm trap (a Rust panic, a read out of bounds, an allocation the heap cannot serve, a stack overflow out of a
 * deep call) may leave the instance's heap broken, so the trapped instance is never called again. Inside the call that
 * trapped:
 *
 * 1. the trap's note is read (crux-core decision wasm-trap-notes: the panic's message and place; with nothing noted,
 *    the glue names the trap itself, an abort without a panic as "ran out of memory (inferred ...)");
 * 2. a fresh instance is made from the compiled module (the glue's `freshInstance()`, crux-core
 *    tools/wasm_glue_fresh.mjs), and the module's settings are made on it again (`settings`);
 * 3. the call is made once more: the objects it uses (the molecule it is a method of, those it is given) are rebuilt
 *    first, each from its journal into the same JS object (`adoptInto`), and every other object at its own first call;
 * 4. if it traps again there, the engine is made fresh once more and the call throws {@link CruxRefused}, a plain
 *    Error saying {@link REFUSED_MESSAGE}: that input alone fails, as a molecule Crux cannot read does, and Crux
 *    stays on.
 *
 * The journal: each object keeps how it was made (a module function, a constructor, or a method of another object at
 * that one's state then, with the arguments) and every call that changed it since (`cruxChanging`: the methods the
 * crate's Rust takes `&mut self` or `self` for, read from its sources at the build; a class it does not list counts as
 * changed by every method). The engine is deterministic, so the journal replayed on a fresh instance makes the same
 * object. Classes in `unjournaled` keep none (a worker's search collections: their molecules are many, and the worker
 * builds them again): their objects are of the old instance after a recovery, and their methods throw (the glue's
 * "...was replaced"), which their owner handles.
 *
 * One console line per trap: the call, the engine, the input's size, the note and what was done; with the debug switch
 * on, a second line with the input itself (structures can be proprietary: off by default). What a user may see (a
 * refused call's error, in a panel or a balloon) names no engine, no call and no input: those are in the console
 * line.
 */
import {isWasmTrap} from './crux-api';

/** How many times a call that trapped is tried again, each time on a fresh instance. */
export const RETRIES = 1;
/** The longest text the debug line quotes in full. */
const DEBUG_TEXT_LIMIT = 100000;
/**
 * What a refused call's error says (a user may see it): no engine, no call, no input. Not a trap's words: a
 * worker's job reports it as a refusal, and CruxService keeps its workers (crux-service.ts `isTrap` reads the message).
 */
export const REFUSED_MESSAGE = 'This structure could not be processed';

/** What a recovery tells: the console line and the listeners ({@link Recovering.onRecover}). */
export interface CruxRecovery {
  /** The call that trapped: `get_mol`, `JSMol.get_svg`, `new ChemReactionScheme`, ... */
  readonly call: string;
  /** The trap as JS saw it: `RuntimeError: unreachable`. */
  readonly trap: string;
  /** The module's note of it (the panic's message and place; else the trap named from itself, an abort without a panic
   * as "ran out of memory (inferred ...)"), or null. */
  readonly note: string | null;
  /** `retried`: answered on a fresh instance; `refused`: it trapped there too, that call failed alone; `recovered`: no
   * call to try again (a free). */
  readonly outcome: 'retried' | 'refused' | 'recovered';
  /** The instance in use now: 0 for the first, one more for each fresh one. */
  readonly instance: number;
  /** The objects rebuilt for the call. */
  readonly rebuilt: number;
  /** From the trap to the answer, ms. */
  readonly ms: number;
}

/**
 * What a call throws when it trapped again on a fresh instance: that input fails; the module is whole. A plain Error
 * saying {@link REFUSED_MESSAGE}, which names no trap (the call and the note are in the console line and in
 * `recovery`).
 */
export class CruxRefused extends Error {
  readonly recovery: CruxRecovery;
  constructor(recovery: CruxRecovery) {
    super(REFUSED_MESSAGE);
    this.name = 'CruxRefused';
    this.recovery = recovery;
  }
}

/** What the patched glue gives (crux-core tools/wasm_glue_fresh.mjs). */
interface FreshGlue {
  freshInstance?(): unknown; // the new instance's exports (InitOutput)
  adoptInto?(target: object, fresh: object): void;
  trapNote?(thrown?: unknown): string | null;
  setTrapHook?(hook: ((thrown: unknown) => void) | null): void;
  instanceNumber?(): number;
  /** Tests: the n-th engine call from now traps in wasm-panic's `cruxTestTrap` (a real panic); 0 disarms. */
  cruxTestArmTrap?(calls: number, panics?: boolean): void;
  cruxChanging?: Readonly<Record<string, readonly string[]>>;
  version?(): string;
}

/** The module, recovering, and what its recovery offers. */
export interface Recovering<T> {
  /** The glue's exports, every call recovering: functions wrapped, classes made recovering (constructors too). */
  crux: T;
  /** Listens to recoveries, told once the call that trapped is over; returns the unsubscribe. */
  onRecover(listener: (recovery: CruxRecovery) => void): () => void;
  /** The trap no fresh instance could get past (the glue has no `freshInstance`, or it trapped twice), or null. */
  failed(): unknown;
  /** Tests: crux-core's real test trap: the `calls`-th engine call from now panics, here and on the next
   * `instances - 1` fresh instances. */
  armTestTrap(calls: number, instances?: number): void;
  /** The debug switch: when on, a trap's second console line gives the call's input; returns it as it is now. */
  debug(on?: boolean): boolean;
  /** How many fresh instances this module has had (0 before the first trap), read at once: a worker compares it
   * around a job to know that its unjournaled objects (the search collections) were lost in it. */
  renewals(): number;
}

export interface RecoverOptions {
  /** What the module sets after its init (`prefer_coordgen(false)`, ...): made again on each fresh instance, given the
   * plain glue and the new instance's exports (its `memory`). */
  settings: (glue: any, exports: {memory?: WebAssembly.Memory}) => void;
  /** How the console lines start: `Chem |`, `Chem | Crux worker |`. */
  log: string;
  /** Classes whose objects keep no journal (see the module's header). */
  unjournaled?: readonly string[];
}

interface Ref {
  readonly journal: Journal;
  readonly upto: number;
}

class HandleArg {
  constructor(readonly ref: Ref) {}
}

interface Op {
  readonly name: string;
  readonly args: readonly unknown[];
}

type Origin =
  | {readonly made: 'function', readonly name: string, readonly args: readonly unknown[]}
  | {readonly made: 'new', readonly cls: Function, readonly args: readonly unknown[]}
  | {readonly made: 'method', readonly name: string, readonly from: Ref, readonly args: readonly unknown[]};

interface Journal {
  readonly seq: number;
  readonly cls: string;
  readonly origin: Origin;
  readonly ops: Op[];
  built: number;
}

type Fn = (...args: unknown[]) => unknown;
type Kind = 'function' | 'static' | 'method';

/** Whether `f` is a class of the glue (a JSMol, a Collection, ...): its prototype has methods. */
const isClass = (f: Function): boolean =>
  !!f.prototype && Object.getOwnPropertyNames(f.prototype).some((name) => name !== 'constructor');

/** The freeing methods: `free`, and chem-wasm's RDKit.js name for it. */
const FREES = new Set(['free', 'delete']);

/**
 * `glue` recovering ({@link Recovering}). Its classes' methods and statics are made recovering in place, once (the
 * objects the module returns are of those classes, wherever they are used); its functions and constructors in the copy
 * of the namespace it returns (`crux`).
 */
export function recoverable<T extends object>(glue: T, options: RecoverOptions): Recovering<T> {
  const fresh = glue as unknown as FreshGlue;
  const classes = new Map<Function, string>();
  for (const [name, value] of Object.entries(glue)) {
    if (typeof value === 'function' && isClass(value))
      classes.set(value, name);
  }
  const changing = new Map<string, Set<string> | null>();
  for (const name of classes.values())
    changing.set(name, fresh.cruxChanging?.[name] ? new Set(fresh.cruxChanging[name]) : null);
  const unjournaled = new Set(options.unjournaled ?? []);

  let dead: unknown = null;
  let pending: {thrown: unknown, in: string} | null = null;
  let debugOn = false;
  let armed = 0;
  let armedInstances = 0;
  let seq = 0;
  let generation = 0;
  let states: States | null = null;
  let rebuiltNow = 0;
  /** The engine's version for the console line, read once: at the first renewal, on the new instance before anything
   * else calls it (a worker wraps its module before instantiating it). Never read in {@link tell}: the glue counts
   * `version()` as an engine call, so a read there, after a renewal, ran on the fresh instance outside any recovering
   * call (a trap there was swallowed) and took the test trap armed on it (the main thread's "traps on a fresh
   * instance too" test then saw Crux read the molecule, the line naming "engine unknown"). */
  let engineVersion: string | null = null;
  const journals = new WeakMap<object, Journal>();
  const live = new Map<number, WeakRef<object>>();
  let pruneAt = 64;
  const listeners = new Set<(recovery: CruxRecovery) => void>();
  /** The plain calls, for the replay: each class's methods by name, the module's functions by name. */
  const methods = new Map<string, Map<string, Fn>>();
  const functions = new Map<string, Fn>();
  const now = (): number => performance.now();

  const classOf = (v: unknown): string | null => {
    if (v === null || typeof v !== 'object')
      return null;
    for (const [cls, name] of classes) {
      if (v instanceof cls)
        return name;
    }
    return null;
  };
  const isHandle = (v: unknown): v is object => classOf(v) !== null;

  const keep = (arg: unknown): unknown => {
    if (isHandle(arg)) {
      const journal = journals.get(arg);
      return journal ? new HandleArg({journal, upto: journal.ops.length}) : arg;
    }
    if (Array.isArray(arg))
      return arg.map(keep);
    if (ArrayBuffer.isView(arg))
      return (arg as unknown as {slice(): unknown}).slice();
    if (arg !== null && typeof arg === 'object') {
      try {
        return structuredClone(arg);
      } catch {
        return arg; // a canvas, a function: kept as it is
      }
    }
    return arg;
  };

  const track = (obj: object, origin: Origin): void => {
    const cls = classOf(obj)!;
    if (unjournaled.has(cls))
      return;
    const journal: Journal = {seq: seq++, cls, origin, ops: [], built: generation};
    journals.set(obj, journal);
    live.set(journal.seq, new WeakRef(obj));
    if (live.size >= pruneAt) {
      for (const [k, ref] of live) {
        if (!ref.deref())
          live.delete(k);
      }
      pruneAt = Math.max(64, live.size * 2);
    }
  };

  /** Whether `name` of class `cls` changes its object: the build's list, or any method of a class it does not list. */
  const changes = (cls: string, name: string): boolean => {
    const set = changing.get(cls);
    return set === null || set === undefined || set.has(name);
  };

  /** Records a call that returned or threw no trap: an object it made, a change it made, an object it freed. */
  const record = (self: unknown, name: string, kind: Kind, args: unknown[], out: unknown, label: string): void => {
    if (isHandle(out) && !journals.has(out)) {
      const from = isHandle(self) ? journals.get(self) : undefined;
      if (kind === 'method' && from)
        track(out, {made: 'method', name, from: {journal: from, upto: from.ops.length}, args: args.map(keep)});
      else if (kind !== 'method' && functions.has(label))
        track(out, {made: 'function', name: label, args: args.map(keep)});
    }
    if (kind !== 'method' || !isHandle(self))
      return;
    const journal = journals.get(self);
    if (!journal)
      return;
    if (FREES.has(name))
      live.delete(journal.seq);
    else if (changes(journal.cls, name))
      journal.ops.push({name, args: args.map(keep)});
  };

  const noteNow = (thrown: unknown): string | null => {
    try {
      return fresh.trapNote?.(thrown) ?? null;
    } catch {
      return null;
    }
  };

  const liveOf = (journal: Journal): object | null => {
    const obj = live.get(journal.seq)?.deref();
    return obj && journals.get(obj) === journal && journal.built !== generation ? obj : null;
  };
  const usableOf = (journal: Journal): object | null => {
    const obj = live.get(journal.seq)?.deref();
    return obj && journals.get(obj) === journal && journal.built === generation ? obj : null;
  };

  /** A fresh instance, the settings made on it again; every object is then rebuilt at its first call. */
  const renew = (cause: unknown): void => {
    if (!fresh.freshInstance || !fresh.adoptInto) {
      dead ??= cause;
      throw dead;
    }
    let last: unknown = cause;
    for (let attempt = 0; attempt < 2; attempt++) {
      try {
        const exports = fresh.freshInstance() as {memory?: WebAssembly.Memory};
        generation++;
        states = new States(methods, functions, isHandle, usableOf, liveOf, (journal, made, into) => {
          fresh.adoptInto!(into, made);
          journal.built = generation;
          rebuiltNow++;
        }, (obj) => methods.get(classOf(obj) ?? '')?.get('free')?.call(obj));
        options.settings(glue, exports);
        if (engineVersion === null) {
          try {
            engineVersion = String(fresh.version?.() ?? 'unknown');
          } catch {
            engineVersion = 'unknown'; // the version is only for the line
          }
        }
        if (armedInstances > 0) {
          armedInstances--;
          fresh.cruxTestArmTrap?.(armed);
        }
        return;
      } catch (thrown) {
        if (!isWasmTrap(thrown))
          throw thrown;
        last = thrown;
      }
    }
    dead ??= last;
    throw dead;
  };

  const ensure = (obj: object): void => {
    const journal = journals.get(obj);
    if (!journal || journal.built === generation || !states)
      return;
    try {
      const made = states.of({journal, upto: journal.ops.length});
      if (journal.built !== generation && made !== obj) {
        fresh.adoptInto!(obj, made);
        journal.built = generation;
        rebuiltNow++;
      }
    } finally {
      states.freeTemps();
    }
  };

  const ensureCall = (self: unknown, args: readonly unknown[]): void => {
    if (generation === 0)
      return;
    if (isHandle(self))
      ensure(self);
    for (const a of args) {
      if (isHandle(a))
        ensure(a);
      else if (Array.isArray(a)) {
        for (const x of a) {
          if (isHandle(x))
            ensure(x);
        }
      }
    }
  };

  const tell = (recovery: CruxRecovery, self: unknown, args: readonly unknown[], thrown: unknown): void => {
    const engine = engineVersion ?? 'unknown';
    const what = recovery.outcome === 'retried' ? 'recovered on a fresh instance, and the call answered there' :
      recovery.outcome === 'refused' ? 'it trapped on a fresh instance too: this call fails, Crux stays on' :
        'recovered on a fresh instance';
    console.warn(`${options.log} Crux trap in ${recovery.call} (engine ${engine}; ${describeInput(self, args)}): ` +
      `${recovery.note ?? recovery.trap}; ${what} ` +
      `(${recovery.rebuilt} objects rebuilt for it, ${recovery.ms.toFixed(1)} ms)`);
    if (debugOn) {
      console.warn(`${options.log} the trapped call's input and the trap's stack (debug switch on)`,
        debugInput(self, args), thrown instanceof Error ? thrown.stack : String(thrown));
    }
    queueMicrotask(() => {
      for (const listener of [...listeners]) {
        try {
          listener(recovery);
        } catch {
          // a listener's failure is its own
        }
      }
    });
  };

  const describeTrap = (thrown: unknown): string => thrown instanceof Error ? `${thrown.name}: ${thrown.message}` :
    String(thrown);
  const instanceNow = (): number => fresh.instanceNumber?.() ?? 0;

  const before = (): void => {
    if (dead)
      throw dead;
    if (pending) {
      const {thrown, in: call} = pending;
      pending = null;
      const t0 = now();
      const note = noteNow(thrown);
      renew(thrown);
      tell({call, trap: describeTrap(thrown), note, outcome: 'recovered', instance: instanceNow(), rebuilt: 0,
        ms: now() - t0}, null, [], thrown);
    }
  };

  /** `plain` recovering: the objects it uses rebuilt first; a trap renews the instance and the call is tried again. */
  const recovering = (plain: Fn, name: string, kind: Kind, label: string,
    call: (self: unknown, args: unknown[]) => unknown): Fn => {
    const freeing = kind === 'method' && FREES.has(name);
    // a function, not an arrow: `this` is the object a method is called on (a method's wrapper stands in for it)
    /* eslint-disable no-invalid-this */
    const wrapped = function(this: unknown, ...args: unknown[]): unknown {
      before();
      rebuiltNow = 0;
      let first: unknown = null;
      let note: string | null = null;
      let trappedAt = 0;
      for (let traps = 0; ;) {
        let called = false;
        try {
          if (!freeing)
            ensureCall(this, args);
          called = true;
          const out = call(this, args);
          record(this, name, kind, args, out, label);
          if (traps) {
            tell({call: label, trap: describeTrap(first), note, outcome: 'retried', instance: instanceNow(),
              rebuilt: rebuiltNow, ms: now() - trappedAt}, this, args, first);
          }
          return out;
        } catch (thrown) {
          if (!isWasmTrap(thrown)) {
            if (called)
              record(this, name, kind, args, undefined, label); // a refusal may have changed part of the object
            if (traps) {
              tell({call: label, trap: describeTrap(first), note, outcome: 'retried', instance: instanceNow(),
                rebuilt: rebuiltNow, ms: now() - trappedAt}, this, args, first);
            }
            throw thrown;
          }
          traps++;
          first ??= thrown;
          note ??= noteNow(thrown);
          trappedAt ||= now();
          if (freeing) {
            const journal = isHandle(this) ? journals.get(this) : undefined;
            if (journal)
              live.delete(journal.seq);
            renew(thrown);
            tell({call: label, trap: describeTrap(first), note, outcome: 'recovered', instance: instanceNow(),
              rebuilt: 0, ms: now() - trappedAt}, null, [], first);
            return undefined;
          }
          renew(thrown);
          if (traps > RETRIES) {
            const recovery: CruxRecovery = {call: label, trap: describeTrap(first), note, outcome: 'refused',
              instance: instanceNow(), rebuilt: rebuiltNow, ms: now() - trappedAt};
            tell(recovery, this, args, first);
            throw new CruxRefused(recovery);
          }
        }
      }
    };
    /* eslint-enable no-invalid-this */
    Object.defineProperty(wrapped, 'length', {value: plain.length});
    Object.defineProperty(wrapped, 'name', {value: plain.name});
    return wrapped;
  };

  // the classes: methods and statics recovering in place, once (every object of the class, wherever it is used)
  const dispose = (Symbol as {dispose?: symbol}).dispose;
  for (const [cls, className] of classes) {
    const own = new Map<string, Fn>();
    methods.set(className, own);
    const prototype = cls.prototype as Record<string | symbol, unknown>;
    const plainFree = prototype.free;
    for (const [target, kind] of [[cls, 'static'], [prototype, 'method']] as const) {
      for (const name of Object.getOwnPropertyNames(target)) {
        if (name === 'constructor' || name === 'prototype' || name === 'length' || name === 'name' ||
          name.startsWith('__'))
          continue;
        const d = Object.getOwnPropertyDescriptor(target, name);
        if (!d || typeof d.value !== 'function' || !d.writable)
          continue;
        const plain = d.value as Fn;
        if (kind === 'method')
          own.set(name, plain);
        else
          functions.set(`${className}.${name}`, plain);
        Object.defineProperty(target, name, {...d, value: recovering(plain, name, kind, `${className}.${name}`,
          (self, args) => plain.apply(self, args))});
      }
    }
    if (dispose && plainFree && prototype[dispose] === plainFree)
      prototype[dispose] = prototype.free;
  }

  // the namespace's copy: functions recovering; each class behind a proxy whose `new` recovers and is journaled (its
  // statics, prototype and `instanceof` are the class's own)
  const copy: Record<string, unknown> = {};
  for (const [name, value] of Object.entries(glue)) {
    if (typeof value !== 'function')
      copy[name] = value;
    else if (classes.has(value)) {
      const base = value as new (...args: unknown[]) => object;
      const make = recovering(base as unknown as Fn, '', 'function', `new ${name}`, (_self, args) => {
        const obj = Reflect.construct(base, args);
        track(obj, {made: 'new', cls: base, args: args.map(keep)});
        return obj;
      });
      copy[name] = new Proxy(base, {construct: (_target, args) => make(...args) as object});
    } else {
      functions.set(name, value as Fn);
      copy[name] = recovering(value as Fn, name, 'function', name, (self, args) => (value as Fn).apply(self, args));
    }
  }
  fresh.setTrapHook?.((thrown) => {
    if (isWasmTrap(thrown))
      pending ??= {thrown, in: 'a free by the browser\'s finalizer'};
  });

  const describeInput = (self: unknown, args: readonly unknown[]): string => {
    const sizeOf = (a: unknown): string => typeof a === 'string' ? `${a.length} chars` :
      typeof a === 'number' || typeof a === 'boolean' || a === null || a === undefined ? String(a) :
        isHandle(a) ? `a ${classOf(a)}` : Array.isArray(a) || ArrayBuffer.isView(a) ?
          `${(a as {length: number}).length} items` : typeof a;
    const on = isHandle(self) ? `on a ${classOf(self)}, ` : '';
    return `${on}input ${args.length ? args.map(sizeOf).join(', ') : 'none'}`;
  };
  const debugInput = (self: unknown, args: readonly unknown[]): unknown => {
    const text = (s: string): string => s.length > DEBUG_TEXT_LIMIT ? `${s.slice(0, DEBUG_TEXT_LIMIT)}...` : s;
    const journalOf = (o: unknown): unknown => {
      const j = isHandle(o) ? journals.get(o) : undefined;
      return j ? {made: j.origin.made, args: j.origin.args.map((a) => typeof a === 'string' ? text(a) : a),
        calls: j.ops.map((op) => op.name)} : null;
    };
    return {args: args.map((a) => typeof a === 'string' ? text(a) : isHandle(a) ? journalOf(a) : a),
      on: journalOf(self)};
  };

  return {
    crux: copy as T,
    onRecover: (listener) => {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    failed: () => dead,
    armTestTrap: (calls, instances = 1) => {
      armed = calls;
      armedInstances = calls > 0 ? Math.max(0, instances - 1) : 0;
      fresh.cruxTestArmTrap?.(calls);
    },
    debug: (on) => {
      if (on !== undefined)
        debugOn = on;
      return debugOn;
    },
    renewals: () => generation,
  };
}

/** One rebuild's states on the instance in use (see crux-sketch's packages/engine/src/recover.ts, the same walk). */
class States {
  readonly #made = new Map<Journal, Map<number, object>>();
  readonly #temps = new Set<object>();

  constructor(
    readonly methods: ReadonlyMap<string, ReadonlyMap<string, Fn>>,
    readonly functions: ReadonlyMap<string, Fn>,
    readonly isHandle: (v: unknown) => v is object,
    readonly usable: (journal: Journal) => object | null,
    readonly pendingOf: (journal: Journal) => object | null,
    readonly adopt: (journal: Journal, made: object, into: object) => void,
    /** Frees a temporary object by its class's plain free. */
    readonly release: (obj: object) => void,
  ) {}

  #hit(ref: Ref): object | undefined {
    const made = this.#made.get(ref.journal)?.get(ref.upto);
    if (made)
      return made;
    return (ref.upto === ref.journal.ops.length && this.usable(ref.journal)) || undefined;
  }

  of(ref: Ref): object {
    const hit = this.#hit(ref);
    if (hit)
      return hit;
    const chain: Ref[] = [];
    let base: object | undefined;
    for (let at: Ref = ref; ;) {
      base = this.#hit(at);
      if (base)
        break;
      chain.push(at);
      if (at.journal.origin.made !== 'method')
        break;
      at = at.journal.origin.from;
    }
    let obj = base;
    for (let i = chain.length - 1; i >= 0; i--) {
      const level = chain[i];
      const origin = level.journal.origin;
      const args = origin.args.map((a) => this.#arg(a));
      const made = origin.made === 'function' ? this.#fn(origin.name)(...args) :
        origin.made === 'new' ? Reflect.construct(origin.cls as new (...a: unknown[]) => object, args) :
          this.#method(level.journal.cls, origin.name, obj!, args, origin.from.journal.cls);
      if (!this.isHandle(made))
        throw new Error(`the recovery's journal: its ${origin.made} made no object`);
      this.#temps.add(made);
      for (const op of level.journal.ops.slice(0, level.upto)) {
        try {
          this.#method(level.journal.cls, op.name, made, op.args.map((a) => this.#arg(a)));
        } catch (thrown) {
          if (isWasmTrap(thrown))
            throw thrown; // a refusal the first time too is replayed as one
        }
      }
      const into = level.upto === level.journal.ops.length ? this.pendingOf(level.journal) : null;
      if (into) {
        this.adopt(level.journal, made, into);
        this.#temps.delete(made);
        obj = into;
      } else {
        let byUpto = this.#made.get(level.journal);
        if (!byUpto)
          this.#made.set(level.journal, byUpto = new Map());
        byUpto.set(level.upto, made);
        obj = made;
      }
    }
    return obj!;
  }

  freeTemps(): void {
    for (const t of this.#temps) {
      try {
        if ((t as {__wbg_ptr?: number}).__wbg_ptr !== 0)
          this.release(t);
      } catch {
        // left to the finalizer
      }
    }
    this.#temps.clear();
    this.#made.clear();
  }

  #arg(a: unknown): unknown {
    if (a instanceof HandleArg)
      return this.of(a.ref);
    if (Array.isArray(a))
      return a.map((x) => this.#arg(x));
    return a;
  }

  #fn(name: string): Fn {
    const fn = this.functions.get(name);
    if (!fn)
      throw new Error(`the recovery's journal names ${name}, which the module does not have`);
    return fn;
  }

  /** A plain method of `cls` (or of `fromCls`, the class of the object a copy was made from) on `self`. */
  #method(cls: string, name: string, self: object, args: unknown[], fromCls?: string): unknown {
    const fn = this.methods.get(fromCls ?? cls)?.get(name);
    if (!fn)
      throw new Error(`the recovery's journal names ${fromCls ?? cls}.${name}, which the module does not have`);
    return fn.apply(self, args);
  }
}
