// An named symbol/brand for detecting Signal instances even when they weren't
// created using the same signals library version.
const BRAND_SYMBOL = Symbol.for("preact-signals");
// Flags for Computed and Effect.
const RUNNING = 1 << 0;
const NOTIFIED = 1 << 1;
const OUTDATED = 1 << 2;
const DISPOSED = 1 << 3;
const HAS_ERROR = 1 << 4;
const TRACKING = 1 << 5;
function startBatch() {
    batchDepth++;
}
function endBatch() {
    if (batchDepth > 1) {
        batchDepth--;
        return;
    }
    let error;
    let hasError = false;
    reconcileBatchSnapshots();
    while (batchedEffect !== undefined) {
        let effect = batchedEffect;
        batchedEffect = undefined;
        batchIteration++;
        while (effect !== undefined) {
            const next = effect._nextBatchedEffect;
            effect._nextBatchedEffect = undefined;
            effect._flags &= ~NOTIFIED;
            if (!(effect._flags & DISPOSED) && needsToRecompute(effect)) {
                try {
                    effect._callback();
                }
                catch (err) {
                    if (!hasError) {
                        error = err;
                        hasError = true;
                    }
                }
            }
            effect = next;
        }
    }
    batchIteration = 0;
    batchDepth--;
    if (hasError) {
        throw error;
    }
}
/**
 * Combine multiple value updates into one "commit" at the end of the provided callback.
 *
 * Batches can be nested and changes are only flushed once the outermost batch callback
 * completes.
 *
 * Accessing a signal that has been modified within a batch will reflect its updated
 * value.
 *
 * @param fn The callback function.
 * @returns The value returned by the callback.
 */
function batch(fn) {
    if (batchDepth > 0) {
        return fn();
    }
    currentBatchSnapshotVersion = ++batchSnapshotVersion;
    /*@__INLINE__**/ startBatch();
    try {
        return fn();
    }
    finally {
        endBatch();
    }
}
// Currently evaluated computed or effect.
let evalContext = undefined;
// Effects captured while constructing a model instance.
let capturedEffects;
/**
 * Run a callback function that can access signal values without
 * subscribing to the signal updates.
 *
 * When called inside a `createModel` factory, this also suppresses
 * model-owned effect capture. Effects created inside the callback will not
 * be owned by the surrounding model and must be disposed manually. Nested
 * `createModel` calls inside the callback still capture their own effects.
 *
 * @param fn The callback function.
 * @returns The value returned by the callback.
 */
function untracked(fn) {
    const prevContext = evalContext;
    const prevCapturedEffects = capturedEffects;
    evalContext = undefined;
    // Model effect capture is another kind of ambient tracking. Suppress it in
    // untracked callbacks while still allowing nested createModel() calls to
    // establish their own capture scope.
    capturedEffects = undefined;
    try {
        return fn();
    }
    finally {
        evalContext = prevContext;
        capturedEffects = prevCapturedEffects;
    }
}
// Effects collected into a batch.
let batchedEffect = undefined;
let batchDepth = 0;
let batchIteration = 0;
let batchSnapshotVersion = 0;
let currentBatchSnapshotVersion = 0;
let batchSnapshots = undefined;
// A global version number for signals, used for fast-pathing repeated
// computed.peek()/computed.value calls when nothing has changed globally.
let globalVersion = 0;
function recordBatchSnapshot(source) {
    // Only capture writes during the user-visible batch callback, not during effect flush.
    if (batchDepth === 0 || batchIteration !== 0) {
        return;
    }
    if (source._batchSnapshotVersion !== currentBatchSnapshotVersion) {
        source._batchSnapshotVersion = currentBatchSnapshotVersion;
        batchSnapshots = {
            _source: source,
            _value: source._value,
            _version: source._version,
            _next: batchSnapshots,
        };
    }
}
function reconcileBatchSnapshots() {
    let snapshots = batchSnapshots;
    batchSnapshots = undefined;
    while (snapshots !== undefined) {
        const source = snapshots._source;
        if (source._value === snapshots._value) {
            // The value was reverted to its pre-batch state. Version numbers must
            // stay monotonic: a lazy computed may have observed an intermediate
            // version during the batch, and rolling the version back would let a
            // future write re-mint that observed number for a different value,
            // making the computed treat it as unchanged forever. Instead,
            // fast-forward subscribers that last saw the pre-batch version so
            // they skip recomputing for the no-op change.
            for (let node = source._targets; node !== undefined; node = node._nextTarget) {
                if (node._version === snapshots._version) {
                    node._version = source._version;
                }
            }
        }
        snapshots = snapshots._next;
    }
}
function addDependency(signal) {
    if (evalContext === undefined) {
        return undefined;
    }
    let node = signal._node;
    if (node === undefined || node._target !== evalContext) {
        /**
         * `signal` is a new dependency. Create a new dependency node, and set it
         * as the tail of the current context's dependency list. e.g:
         *
         * { A <-> B       }
         *         ↑     ↑
         *        tail  node (new)
         *               ↓
         * { A <-> B <-> C }
         *               ↑
         *              tail (evalContext._sources)
         */
        node = {
            _version: 0,
            _source: signal,
            _prevSource: evalContext._sources,
            _nextSource: undefined,
            _target: evalContext,
            _prevTarget: undefined,
            _nextTarget: undefined,
            _rollbackNode: node,
        };
        if (evalContext._sources !== undefined) {
            evalContext._sources._nextSource = node;
        }
        evalContext._sources = node;
        signal._node = node;
        // Subscribe to change notifications from this dependency if we're in an effect
        // OR evaluating a computed signal that in turn has subscribers.
        if (evalContext._flags & TRACKING) {
            signal._subscribe(node);
        }
        return node;
    }
    else if (node._version === -1) {
        // `signal` is an existing dependency from a previous evaluation. Reuse it.
        node._version = 0;
        /**
         * If `node` is not already the current tail of the dependency list (i.e.
         * there is a next node in the list), then make the `node` the new tail. e.g:
         *
         * { A <-> B <-> C <-> D }
         *         ↑           ↑
         *        node   ┌─── tail (evalContext._sources)
         *         └─────│─────┐
         *               ↓     ↓
         * { A <-> C <-> D <-> B }
         *                     ↑
         *                    tail (evalContext._sources)
         */
        if (node._nextSource !== undefined) {
            node._nextSource._prevSource = node._prevSource;
            if (node._prevSource !== undefined) {
                node._prevSource._nextSource = node._nextSource;
            }
            node._prevSource = evalContext._sources;
            node._nextSource = undefined;
            evalContext._sources._nextSource = node;
            evalContext._sources = node;
        }
        // We can assume that the currently evaluated effect / computed signal is already
        // subscribed to change notifications from `signal` if needed.
        return node;
    }
    return undefined;
}
/** @internal */
// A class with the same name has already been declared, so we need to ignore
// TypeScript's warning about a redeclared variable.
//
// The previously declared class is implemented here with ES5-style prototypes.
// This enables better control of the transpiled output size.
// @ts-ignore: "Cannot redeclare exported variable 'Signal'."
function Signal(value, options) {
    this._value = value;
    this._version = 0;
    this._node = undefined;
    this._targets = undefined;
    this._batchSnapshotVersion = 0;
    this._watched = options === null || options === void 0 ? void 0 : options.watched;
    this._unwatched = options === null || options === void 0 ? void 0 : options.unwatched;
    this.name = options === null || options === void 0 ? void 0 : options.name;
}
Signal.prototype.brand = BRAND_SYMBOL;
Signal.prototype._refresh = function () {
    return true;
};
Signal.prototype._subscribe = function (node) {
    const targets = this._targets;
    if (targets !== node && node._prevTarget === undefined) {
        node._nextTarget = targets;
        this._targets = node;
        if (targets !== undefined) {
            targets._prevTarget = node;
        }
        else {
            untracked(() => {
                var _a;
                (_a = this._watched) === null || _a === void 0 ? void 0 : _a.call(this);
            });
        }
    }
};
Signal.prototype._unsubscribe = function (node) {
    // Only run the unsubscribe step if the signal has any subscribers to begin with.
    if (this._targets !== undefined) {
        const prev = node._prevTarget;
        const next = node._nextTarget;
        if (prev !== undefined) {
            prev._nextTarget = next;
            node._prevTarget = undefined;
        }
        if (next !== undefined) {
            next._prevTarget = prev;
            node._nextTarget = undefined;
        }
        if (node === this._targets) {
            this._targets = next;
            if (next === undefined) {
                untracked(() => {
                    var _a;
                    (_a = this._unwatched) === null || _a === void 0 ? void 0 : _a.call(this);
                });
            }
        }
    }
};
Signal.prototype.subscribe = function (fn) {
    return effect(() => {
        const value = this.value;
        untracked(() => fn(value));
    }, { name: "sub" });
};
Signal.prototype.valueOf = function () {
    return this.value;
};
Signal.prototype.toString = function () {
    return this.value + "";
};
Signal.prototype.toJSON = function () {
    return this.value;
};
Signal.prototype.peek = function () {
    return untracked(() => this.value);
};
Object.defineProperty(Signal.prototype, "value", {
    get() {
        const node = addDependency(this);
        if (node !== undefined) {
            node._version = this._version;
        }
        return this._value;
    },
    set(value) {
        if (value !== this._value) {
            if (batchIteration > 100) {
                throw new Error("Cycle detected");
            }
            recordBatchSnapshot(this);
            this._value = value;
            this._version++;
            globalVersion++;
            /**@__INLINE__*/ startBatch();
            try {
                for (let node = this._targets; node !== undefined; node = node._nextTarget) {
                    node._target._notify();
                }
            }
            finally {
                endBatch();
            }
        }
    },
});
export function signal(value, options) {
    return new Signal(value, options);
}
//#endregion Signal
//#region Computed
function needsToRecompute(target) {
    // Check the dependencies for changed values. The dependency list is already
    // in order of use. Therefore if multiple dependencies have changed values, only
    // the first used dependency is re-evaluated at this point.
    for (let node = target._sources; node !== undefined; node = node._nextSource) {
        if (
        // If the dependency has definitely been updated since its version number
        // was observed, then we need to recompute. This first check is not strictly
        // necessary for correctness, but allows us to skip the refresh call if the
        // dependency has already been updated.
        node._source._version !== node._version ||
            // Refresh the dependency. If there's something blocking the refresh (e.g. a
            // dependency cycle), then we need to recompute.
            !node._source._refresh() ||
            // If the dependency got a new version after the refresh, then we need to recompute.
            node._source._version !== node._version) {
            return true;
        }
    }
    // If none of the dependencies have changed values since last recompute then
    // there's no need to recompute.
    return false;
}
function prepareSources(target) {
    /**
     * 1. Mark all current sources as re-usable nodes (version: -1)
     * 2. Set a rollback node if the current node is being used in a different context
     * 3. Point 'target._sources' to the tail of the doubly-linked list, e.g:
     *
     *    { undefined <- A <-> B <-> C -> undefined }
     *                   ↑           ↑
     *                   │           └──────┐
     * target._sources = A; (node is head)  │
     *                   ↓                  │
     * target._sources = C; (node is tail) ─┘
     */
    for (let node = target._sources; node !== undefined; node = node._nextSource) {
        const rollbackNode = node._source._node;
        if (rollbackNode !== undefined) {
            node._rollbackNode = rollbackNode;
        }
        node._source._node = node;
        node._version = -1;
        if (node._nextSource === undefined) {
            target._sources = node;
            break;
        }
    }
}
function cleanupSources(target) {
    let node = target._sources;
    let head = undefined;
    /**
     * At this point 'target._sources' points to the tail of the doubly-linked list.
     * It contains all existing sources + new sources in order of use.
     * Iterate backwards until we find the head node while dropping old dependencies.
     */
    while (node !== undefined) {
        const prev = node._prevSource;
        /**
         * The node was not re-used, unsubscribe from its change notifications and remove itself
         * from the doubly-linked list. e.g:
         *
         * { A <-> B <-> C }
         *         ↓
         *    { A <-> C }
         */
        if (node._version === -1) {
            node._source._unsubscribe(node);
            if (prev !== undefined) {
                prev._nextSource = node._nextSource;
            }
            if (node._nextSource !== undefined) {
                node._nextSource._prevSource = prev;
            }
        }
        else {
            /**
             * The new head is the last node seen which wasn't removed/unsubscribed
             * from the doubly-linked list. e.g:
             *
             * { A <-> B <-> C }
             *   ↑     ↑     ↑
             *   │     │     └ head = node
             *   │     └ head = node
             *   └ head = node
             */
            head = node;
        }
        node._source._node = node._rollbackNode;
        if (node._rollbackNode !== undefined) {
            node._rollbackNode = undefined;
        }
        node = prev;
    }
    target._sources = head;
}
/** @internal */
function Computed(fn, options) {
    Signal.call(this, undefined, options);
    this._fn = fn;
    this._sources = undefined;
    this._globalVersion = globalVersion - 1;
    this._flags = OUTDATED;
}
Computed.prototype = new Signal();
Computed.prototype._refresh = function () {
    this._flags &= ~NOTIFIED;
    if (this._flags & RUNNING) {
        return false;
    }
    // If this computed signal has subscribed to updates from its dependencies
    // (TRACKING flag set) and none of them have notified about changes (OUTDATED
    // flag not set), then the computed value can't have changed.
    if ((this._flags & (OUTDATED | TRACKING)) === TRACKING) {
        return true;
    }
    this._flags &= ~OUTDATED;
    if (this._globalVersion === globalVersion) {
        return true;
    }
    this._globalVersion = globalVersion;
    // Mark this computed signal running before checking the dependencies for value
    // changes, so that the RUNNING flag can be used to notice cyclical dependencies.
    this._flags |= RUNNING;
    if (this._version > 0 && !needsToRecompute(this)) {
        this._flags &= ~RUNNING;
        return true;
    }
    const prevContext = evalContext;
    try {
        prepareSources(this);
        evalContext = this;
        const value = this._fn();
        if (this._flags & HAS_ERROR ||
            this._value !== value ||
            this._version === 0) {
            this._value = value;
            this._flags &= ~HAS_ERROR;
            this._version++;
        }
    }
    catch (err) {
        this._value = err;
        this._flags |= HAS_ERROR;
        this._version++;
    }
    evalContext = prevContext;
    cleanupSources(this);
    this._flags &= ~RUNNING;
    return true;
};
Computed.prototype._subscribe = function (node) {
    if (this._targets === undefined) {
        this._flags |= OUTDATED | TRACKING;
        // A computed signal subscribes lazily to its dependencies when it
        // gets its first subscriber.
        for (let node = this._sources; node !== undefined; node = node._nextSource) {
            node._source._subscribe(node);
        }
    }
    Signal.prototype._subscribe.call(this, node);
};
Computed.prototype._unsubscribe = function (node) {
    // Only run the unsubscribe step if the computed signal has any subscribers.
    if (this._targets !== undefined) {
        Signal.prototype._unsubscribe.call(this, node);
        // Computed signal unsubscribes from its dependencies when it loses its last subscriber.
        // This makes it possible for unreferences subgraphs of computed signals to get garbage collected.
        if (this._targets === undefined) {
            this._flags &= ~TRACKING;
            for (let node = this._sources; node !== undefined; node = node._nextSource) {
                node._source._unsubscribe(node);
            }
        }
    }
};
Computed.prototype._notify = function () {
    if (!(this._flags & NOTIFIED)) {
        this._flags |= OUTDATED | NOTIFIED;
        for (let node = this._targets; node !== undefined; node = node._nextTarget) {
            node._target._notify();
        }
    }
};
Object.defineProperty(Computed.prototype, "value", {
    get() {
        if (this._flags & RUNNING) {
            throw new Error("Cycle detected");
        }
        const node = addDependency(this);
        this._refresh();
        if (node !== undefined) {
            node._version = this._version;
        }
        if (this._flags & HAS_ERROR) {
            throw this._value;
        }
        return this._value;
    },
});
/**
 * Create a new signal that is computed based on the values of other signals.
 *
 * The returned computed signal is read-only, and its value is automatically
 * updated when any signals accessed from within the callback function change.
 *
 * @param fn The effect callback.
 * @returns A new read-only signal.
 */
function computed(fn, options) {
    return new Computed(fn, options);
}
//#endregion Computed
//#region Effect
function cleanupEffect(effect) {
    const cleanup = effect._cleanup;
    effect._cleanup = undefined;
    if (typeof cleanup === "function") {
        /*@__INLINE__**/ startBatch();
        // Run cleanup functions always outside of any context.
        const prevContext = evalContext;
        evalContext = undefined;
        try {
            cleanup();
        }
        catch (err) {
            effect._flags &= ~RUNNING;
            effect._flags |= DISPOSED;
            disposeEffect(effect);
            throw err;
        }
        finally {
            evalContext = prevContext;
            endBatch();
        }
    }
}
function disposeEffect(effect) {
    for (let node = effect._sources; node !== undefined; node = node._nextSource) {
        node._source._unsubscribe(node);
    }
    effect._fn = undefined;
    effect._sources = undefined;
    cleanupEffect(effect);
}
function endEffect(prevContext) {
    if (evalContext !== this) {
        throw new Error("Out-of-order effect");
    }
    cleanupSources(this);
    evalContext = prevContext;
    this._flags &= ~RUNNING;
    if (this._flags & DISPOSED) {
        disposeEffect(this);
    }
    endBatch();
}
/** @internal */
function Effect(fn, options) {
    this._fn = fn;
    this._cleanup = undefined;
    this._sources = undefined;
    this._nextBatchedEffect = undefined;
    this._flags = TRACKING;
    this.name = options === null || options === void 0 ? void 0 : options.name;
    if (capturedEffects) {
        capturedEffects.push(this);
    }
}
Effect.prototype._callback = function () {
    const finish = this._start();
    try {
        if (this._flags & DISPOSED)
            return;
        if (this._fn === undefined)
            return;
        const cleanup = this._fn();
        if (typeof cleanup === "function") {
            this._cleanup = cleanup;
        }
    }
    finally {
        finish();
    }
};
Effect.prototype._start = function () {
    if (this._flags & RUNNING) {
        throw new Error("Cycle detected");
    }
    this._flags |= RUNNING;
    this._flags &= ~DISPOSED;
    cleanupEffect(this);
    prepareSources(this);
    /*@__INLINE__**/ startBatch();
    const prevContext = evalContext;
    evalContext = this;
    return endEffect.bind(this, prevContext);
};
Effect.prototype._notify = function () {
    if (!(this._flags & NOTIFIED)) {
        this._flags |= NOTIFIED;
        this._nextBatchedEffect = batchedEffect;
        batchedEffect = this;
    }
};
Effect.prototype._dispose = function () {
    this._flags |= DISPOSED;
    if (!(this._flags & RUNNING)) {
        disposeEffect(this);
    }
};
Effect.prototype.dispose = function () {
    this._dispose();
};
/**
 * Create an effect to run arbitrary code in response to signal changes.
 *
 * An effect tracks which signals are accessed within the given callback
 * function `fn`, and re-runs the callback when those signals change.
 *
 * The callback may return a cleanup function. The cleanup function gets
 * run once, either when the callback is next called or when the effect
 * gets disposed, whichever happens first.
 *
 * @param fn The effect callback.
 * @returns A function for disposing the effect.
 */
function effect(fn, options) {
    const effect = new Effect(fn, options);
    try {
        effect._callback();
    }
    catch (err) {
        effect._dispose();
        throw err;
    }
    // Return a bound function instead of a wrapper like `() => effect._dispose()`,
    // because bound functions seem to be just as fast and take up a lot less memory.
    const dispose = effect._dispose.bind(effect);
    dispose[Symbol.dispose] = dispose;
    return dispose;
}
//#endregion Effect
//#region Action
function action(fn) {
    return function actionWrapper(...args) {
        return batch(() => untracked(() => fn.apply(this, args)));
    };
}
function startCapturingEffects() {
    let prevCapturedEffects = capturedEffects;
    // Always establish a fresh capture scope, even when `untracked()` has
    // temporarily cleared the parent scope. This lets nested models own their
    // effects without promoting them to a suppressed outer scope.
    capturedEffects = [];
    return function stopCapturingEffects() {
        let modelEffects = capturedEffects;
        if (capturedEffects && prevCapturedEffects) {
            prevCapturedEffects = prevCapturedEffects.concat(capturedEffects);
        }
        capturedEffects = prevCapturedEffects;
        return modelEffects;
    };
}
const wrapInAction = (value) => {
    for (const key in value) {
        const val = value[key];
        if (typeof val === "function") {
            value[key] = action(val);
        }
        else if (typeof val === "object" && val !== null && !("brand" in val)) {
            // Recursively wrap nested object properties in actions. This allows users to write
            // nested models without worrying about wrapping their functions in `action`.
            wrapInAction(val);
        }
    }
};
function createModel(modelFactory) {
    return function SignalModel(...args) {
        let modelEffects;
        let model;
        const stopCapturingEffects = startCapturingEffects();
        try {
            model = modelFactory(...args);
        }
        catch (err) {
            // Drop any captured effects on error. Errors from nested models will bubble
            // up here and recursively reset `capturedEffects` to `undefined` preventing
            // any captured effects from leaking
            capturedEffects = undefined;
            throw err;
        }
        finally {
            modelEffects = stopCapturingEffects();
        }
        wrapInAction(model);
        model[Symbol.dispose] = action(function disposeModel() {
            if (modelEffects) {
                for (let i = 0; i < modelEffects.length; i++) {
                    modelEffects[i].dispose();
                }
            }
            modelEffects = undefined;
        });
        return model;
    };
}
//#endregion createModel
export { computed, effect, batch, untracked, action, createModel, Signal, Effect, Computed, };
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoic2lnbmFscy1jb3JlLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsic2lnbmFscy1jb3JlLnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBLDhFQUE4RTtBQUM5RSxrREFBa0Q7QUFDbEQsTUFBTSxZQUFZLEdBQUcsTUFBTSxDQUFDLEdBQUcsQ0FBQyxnQkFBZ0IsQ0FBQyxDQUFDO0FBRWxELGlDQUFpQztBQUNqQyxNQUFNLE9BQU8sR0FBRyxDQUFDLElBQUksQ0FBQyxDQUFDO0FBQ3ZCLE1BQU0sUUFBUSxHQUFHLENBQUMsSUFBSSxDQUFDLENBQUM7QUFDeEIsTUFBTSxRQUFRLEdBQUcsQ0FBQyxJQUFJLENBQUMsQ0FBQztBQUN4QixNQUFNLFFBQVEsR0FBRyxDQUFDLElBQUksQ0FBQyxDQUFDO0FBQ3hCLE1BQU0sU0FBUyxHQUFHLENBQUMsSUFBSSxDQUFDLENBQUM7QUFDekIsTUFBTSxRQUFRLEdBQUcsQ0FBQyxJQUFJLENBQUMsQ0FBQztBQTBCeEIsU0FBUyxVQUFVO0lBQ2xCLFVBQVUsRUFBRSxDQUFDO0FBQ2QsQ0FBQztBQUVELFNBQVMsUUFBUTtJQUNoQixJQUFJLFVBQVUsR0FBRyxDQUFDLEVBQUUsQ0FBQztRQUNwQixVQUFVLEVBQUUsQ0FBQztRQUNiLE9BQU87SUFDUixDQUFDO0lBRUQsSUFBSSxLQUFjLENBQUM7SUFDbkIsSUFBSSxRQUFRLEdBQUcsS0FBSyxDQUFDO0lBQ3JCLHVCQUF1QixFQUFFLENBQUM7SUFFMUIsT0FBTyxhQUFhLEtBQUssU0FBUyxFQUFFLENBQUM7UUFDcEMsSUFBSSxNQUFNLEdBQXVCLGFBQWEsQ0FBQztRQUMvQyxhQUFhLEdBQUcsU0FBUyxDQUFDO1FBRTFCLGNBQWMsRUFBRSxDQUFDO1FBRWpCLE9BQU8sTUFBTSxLQUFLLFNBQVMsRUFBRSxDQUFDO1lBQzdCLE1BQU0sSUFBSSxHQUF1QixNQUFNLENBQUMsa0JBQWtCLENBQUM7WUFDM0QsTUFBTSxDQUFDLGtCQUFrQixHQUFHLFNBQVMsQ0FBQztZQUN0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UUFDdkMsQ0FBQztRQUVELFlBQVksQ0FBQyxLQUFLLENBQUMsQ0FBQztRQUVwQixLQUFLLENBQUMsTUFBTSxDQUFDLE9BQU8sQ0FBQyxHQUFHLE1BQU0sQ0FBQyxTQUFTLFlBQVk7WUFDbkQsSUFBSSxZQUFZLEVBQUUsQ0FBQztnQkFDbEIsS0FBSyxJQUFJLENBQUMsR0FBRyxDQUFDLEVBQUUsQ0FBQyxHQUFHLFlBQVksQ0FBQyxNQUFNLEVBQUUsQ0FBQyxFQUFFLEVBQUUsQ0FBQztvQkFDOUMsWUFBWSxDQUFDLENBQUMsQ0FBQyxDQUFDLE9BQU8sRUFBRSxDQUFDO2dCQUMzQixDQUFDO1lBQ0YsQ0FBQztZQUVELFlBQVksR0FBRyxTQUFTLENBQUM7UUFDMUIsQ0FBQyxDQUFDLENBQUM7UUFFSCxPQUFPLEtBQUssQ0FBQztJQUNkLENBQW1ELENBQUM7QUFDckQsQ0FBQztBQUVELHdCQUF3QjtBQUV4QixPQUFPLEVBQ04sUUFBUSxFQUNSLE1BQU0sRUFDTixLQUFLLEVBQ0wsU0FBUyxFQUNULE1BQU0sRUFDTixXQUFXLEVBQ1gsTUFBTSxFQUVOLE1BQU0sRUFDTixRQUFRLEdBQ1IsQ0FBQyIsInNvdXJjZXNDb250ZW50IjpbIi8vIEFuIG5hbWVkIHN5bWJvbC9icmFuZCBmb3IgZGV0ZWN0aW5nIFNpZ25hbCBpbnN0YW5jZXMgZXZlbiB3aGVuIHRoZXkgd2VyZW4ndFxuLy8gY3JlYXRlZCB1c2luZyB0aGUgc2FtZSBzaWduYWxzIGxpYnJhcnkgdmVyc2lvbi5cbmNvbnN0IEJSQU5EX1NZTUJPTCA9IFN5bWJvbC5mb3IoXCJwcmVhY3Qtc2lnbmFsc1wiKTtcblxuLy8gRmxhZ3MgZm9yIENvbXB1dGVkIGFuZCBFZmZlY3QuXG5jb25zdCBSVU5OSU5HID0gMSA8PCAwO1xuY29uc3QgTk9USUZJRUQgPSAxIDw8IDE7XG5jb25zdCBPVVREQVRFRCA9IDEgPDwgMjtcbmNvbnN0IERJU1BPU0VEID0gMSA8PCAzO1xuY29uc3QgSEFTX0VSUk9SID0gMSA8PCA0O1xuY29uc3QgVFJBQ0tJTkcgPSAxIDw8IDU7XG5cbi8vIEEgbGlua2VkIGxpc3Qgbm9kZSB1c2VkIHRvIHRyYWNrIGRlcGVuZGVuY2llcyAoc291cmNlcykgYW5kIGRlcGVuZGVudHMgKHRhcmdldHMpLlxuLy8gQWxzbyB1c2VkIHRvIHJlbWVtYmVyIHRoZSBzb3VyY2UncyBsYXN0IHZlcnNpb24gbnVtYmVyIHRoYXQgdGhlIHRhcmdldCBzYXcuXG50eXBlIE5vZGUgPSB7XG5cdC8vIEEgc291cmNlIHdob3NlIHZhbHVlIHRoZSB0YXJnZXQgZGVwZW5kcyBvbi5cblx0X3NvdXJjZTogU2lnbmFsO1xuXHRfcHJldlNvdXJjZT86IE5vZGU7XG5cdF9uZXh0U291cmNlPzogTm9kZTtcblxuXHQvLyBBIHRhcmdldCB0aGF0IGRlcGVuZHMgb24gdGhlIHNvdXJjZSBhbmQgc2hvdWxkIGJlIG5vdGlmaWVkIHdoZW4gdGhlIHNvdXJjZSBjaGFuZ2VzLlxuXHRfdGFyZ2V0OiBDb21wdXRlZCB8IEVmZmVjdDtcblx0X3ByZXZUYXJnZXQ/OiBOb2RlO1xuXHRfbmV4dFRhcmdldD86IE5vZGU7XG5cblx0Ly8gVGhlIHZlcnNpb24gbnVtYmVyIG9mIHRoZSBzb3VyY2UgdGhhdCB0YXJnZXQgaGFzIGxhc3Qgc2Vlbi4gV2UgdXNlIHZlcnNpb24gbnVtYmVyc1xuXHQvLyBpbnN0ZWFkIG9mIHN0b3JpbmcgdGhlIHNvdXJjZSB2YWx1ZSwgYmVjYXVzZSBzb3VyY2UgdmFsdWVzIGNhbiB0YWtlIGFyYml0cmFyeSBhbW91bnRcblx0Ly8gb2YgbWVtb3J5LCBhbmQgY29tcHV0ZWRzIGNvdWxkIGhhbmcgb24gdG8gdGhlbSBmb3JldmVyIGJlY2F1c2UgdGhleSdyZSBsYXppbHkgZXZhbHVhdGVkLlxuXHQvLyBVc2UgdGhlIHNwZWNpYWwgdmFsdWUgLTEgdG8gbWFyayBwb3RlbnRpYWxseSB1bnVzZWQgYnV0IHJlY3ljbGFibGUgbm9kZXMuXG5cdF92ZXJzaW9uOiBudW1iZXI7XG5cblx0Ly8gVXNlZCB0byByZW1lbWJlciAmIHJvbGwgYmFjayB0aGUgc291cmNlJ3MgcHJldmlvdXMgYC5fbm9kZWAgdmFsdWUgd2hlbiBlbnRlcmluZyAmXG5cdC8vIGV4aXRpbmcgYSBuZXcgZXZhbHVhdGlvbiBjb250ZXh0LlxuXHRfcm9sbGJhY2tOb2RlPzogTm9kZTtcbn07XG5cbmZ1bmN0aW9uIHN0YXJ0QmF0Y2goKSB7XG5cdGJhdGNoRGVwdGgrKztcbn1cblxuZnVuY3Rpb24gZW5kQmF0Y2goKSB7XG5cdGlmIChiYXRjaERlcHRoID4gMSkge1xuXHRcdGJhdGNoRGVwdGgtLTtcblx0XHRyZXR1cm47XG5cdH1cblxuXHRsZXQgZXJyb3I6IHVua25vd247XG5cdGxldCBoYXNFcnJvciA9IGZhbHNlO1xuXHRyZWNvbmNpbGVCYXRjaFNuYXBzaG90cygpO1xuXG5cdHdoaWxlIChiYXRjaGVkRWZmZWN0ICE9PSB1bmRlZmluZWQpIHtcblx0XHRsZXQgZWZmZWN0OiBFZmZlY3QgfCB1bmRlZmluZWQgPSBiYXRjaGVkRWZmZWN0O1xuXHRcdGJhdGNoZWRFZmZlY3QgPSB1bmRlZmluZWQ7XG5cblx0XHRiYXRjaEl0ZXJhdGlvbisrO1xuXG5cdFx0d2hpbGUgKGVmZmVjdCAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRjb25zdCBuZXh0OiBFZmZlY3QgfCB1bmRlZmluZWQgPSBlZmZlY3QuX25leHRCYXRjaGVkRWZmZWN0O1xuXHRcdFx0ZWZmZWN0Ll9uZXh0QmF0Y2hlZEVmZmVjdCA9IHVuZGVmaW5lZDtcblx0XHRcdGVmZmVjdC5fZmxhZ3MgJj0gfk5PVElGSUVEO1xuXG5cdFx0XHRpZiAoIShlZmZlY3QuX2ZsYWdzICYgRElTUE9TRUQpICYmIG5lZWRzVG9SZWNvbXB1dGUoZWZmZWN0KSkge1xuXHRcdFx0XHR0cnkge1xuXHRcdFx0XHRcdGVmZmVjdC5fY2FsbGJhY2soKTtcblx0XHRcdFx0fSBjYXRjaCAoZXJyKSB7XG5cdFx0XHRcdFx0aWYgKCFoYXNFcnJvcikge1xuXHRcdFx0XHRcdFx0ZXJyb3IgPSBlcnI7XG5cdFx0XHRcdFx0XHRoYXNFcnJvciA9IHRydWU7XG5cdFx0XHRcdFx0fVxuXHRcdFx0XHR9XG5cdFx0XHR9XG5cdFx0XHRlZmZlY3QgPSBuZXh0O1xuXHRcdH1cblx0fVxuXHRiYXRjaEl0ZXJhdGlvbiA9IDA7XG5cdGJhdGNoRGVwdGgtLTtcblxuXHRpZiAoaGFzRXJyb3IpIHtcblx0XHR0aHJvdyBlcnJvcjtcblx0fVxufVxuXG4vKipcbiAqIENvbWJpbmUgbXVsdGlwbGUgdmFsdWUgdXBkYXRlcyBpbnRvIG9uZSBcImNvbW1pdFwiIGF0IHRoZSBlbmQgb2YgdGhlIHByb3ZpZGVkIGNhbGxiYWNrLlxuICpcbiAqIEJhdGNoZXMgY2FuIGJlIG5lc3RlZCBhbmQgY2hhbmdlcyBhcmUgb25seSBmbHVzaGVkIG9uY2UgdGhlIG91dGVybW9zdCBiYXRjaCBjYWxsYmFja1xuICogY29tcGxldGVzLlxuICpcbiAqIEFjY2Vzc2luZyBhIHNpZ25hbCB0aGF0IGhhcyBiZWVuIG1vZGlmaWVkIHdpdGhpbiBhIGJhdGNoIHdpbGwgcmVmbGVjdCBpdHMgdXBkYXRlZFxuICogdmFsdWUuXG4gKlxuICogQHBhcmFtIGZuIFRoZSBjYWxsYmFjayBmdW5jdGlvbi5cbiAqIEByZXR1cm5zIFRoZSB2YWx1ZSByZXR1cm5lZCBieSB0aGUgY2FsbGJhY2suXG4gKi9cbmZ1bmN0aW9uIGJhdGNoPFQ+KGZuOiAoKSA9PiBUKTogVCB7XG5cdGlmIChiYXRjaERlcHRoID4gMCkge1xuXHRcdHJldHVybiBmbigpO1xuXHR9XG5cdGN1cnJlbnRCYXRjaFNuYXBzaG90VmVyc2lvbiA9ICsrYmF0Y2hTbmFwc2hvdFZlcnNpb247XG5cdC8qQF9fSU5MSU5FX18qKi8gc3RhcnRCYXRjaCgpO1xuXHR0cnkge1xuXHRcdHJldHVybiBmbigpO1xuXHR9IGZpbmFsbHkge1xuXHRcdGVuZEJhdGNoKCk7XG5cdH1cbn1cblxuLy8gQ3VycmVudGx5IGV2YWx1YXRlZCBjb21wdXRlZCBvciBlZmZlY3QuXG5sZXQgZXZhbENvbnRleHQ6IENvbXB1dGVkIHwgRWZmZWN0IHwgdW5kZWZpbmVkID0gdW5kZWZpbmVkO1xuXG4vLyBFZmZlY3RzIGNhcHR1cmVkIHdoaWxlIGNvbnN0cnVjdGluZyBhIG1vZGVsIGluc3RhbmNlLlxubGV0IGNhcHR1cmVkRWZmZWN0czogRWZmZWN0W10gfCB1bmRlZmluZWQ7XG5cbi8qKlxuICogUnVuIGEgY2FsbGJhY2sgZnVuY3Rpb24gdGhhdCBjYW4gYWNjZXNzIHNpZ25hbCB2YWx1ZXMgd2l0aG91dFxuICogc3Vic2NyaWJpbmcgdG8gdGhlIHNpZ25hbCB1cGRhdGVzLlxuICpcbiAqIFdoZW4gY2FsbGVkIGluc2lkZSBhIGBjcmVhdGVNb2RlbGAgZmFjdG9yeSwgdGhpcyBhbHNvIHN1cHByZXNzZXNcbiAqIG1vZGVsLW93bmVkIGVmZmVjdCBjYXB0dXJlLiBFZmZlY3RzIGNyZWF0ZWQgaW5zaWRlIHRoZSBjYWxsYmFjayB3aWxsIG5vdFxuICogYmUgb3duZWQgYnkgdGhlIHN1cnJvdW5kaW5nIG1vZGVsIGFuZCBtdXN0IGJlIGRpc3Bvc2VkIG1hbnVhbGx5LiBOZXN0ZWRcbiAqIGBjcmVhdGVNb2RlbGAgY2FsbHMgaW5zaWRlIHRoZSBjYWxsYmFjayBzdGlsbCBjYXB0dXJlIHRoZWlyIG93biBlZmZlY3RzLlxuICpcbiAqIEBwYXJhbSBmbiBUaGUgY2FsbGJhY2sgZnVuY3Rpb24uXG4gKiBAcmV0dXJucyBUaGUgdmFsdWUgcmV0dXJuZWQgYnkgdGhlIGNhbGxiYWNrLlxuICovXG5mdW5jdGlvbiB1bnRyYWNrZWQ8VD4oZm46ICgpID0+IFQpOiBUIHtcblx0Y29uc3QgcHJldkNvbnRleHQgPSBldmFsQ29udGV4dDtcblx0Y29uc3QgcHJldkNhcHR1cmVkRWZmZWN0cyA9IGNhcHR1cmVkRWZmZWN0cztcblxuXHRldmFsQ29udGV4dCA9IHVuZGVmaW5lZDtcblx0Ly8gTW9kZWwgZWZmZWN0IGNhcHR1cmUgaXMgYW5vdGhlciBraW5kIG9mIGFtYmllbnQgdHJhY2tpbmcuIFN1cHByZXNzIGl0IGluXG5cdC8vIHVudHJhY2tlZCBjYWxsYmFja3Mgd2hpbGUgc3RpbGwgYWxsb3dpbmcgbmVzdGVkIGNyZWF0ZU1vZGVsKCkgY2FsbHMgdG9cblx0Ly8gZXN0YWJsaXNoIHRoZWlyIG93biBjYXB0dXJlIHNjb3BlLlxuXHRjYXB0dXJlZEVmZmVjdHMgPSB1bmRlZmluZWQ7XG5cdHRyeSB7XG5cdFx0cmV0dXJuIGZuKCk7XG5cdH0gZmluYWxseSB7XG5cdFx0ZXZhbENvbnRleHQgPSBwcmV2Q29udGV4dDtcblx0XHRjYXB0dXJlZEVmZmVjdHMgPSBwcmV2Q2FwdHVyZWRFZmZlY3RzO1xuXHR9XG59XG5cbi8vIEVmZmVjdHMgY29sbGVjdGVkIGludG8gYSBiYXRjaC5cbmxldCBiYXRjaGVkRWZmZWN0OiBFZmZlY3QgfCB1bmRlZmluZWQgPSB1bmRlZmluZWQ7XG5sZXQgYmF0Y2hEZXB0aCA9IDA7XG5sZXQgYmF0Y2hJdGVyYXRpb24gPSAwO1xuXG50eXBlIEJhdGNoU25hcHNob3QgPSB7XG5cdF9zb3VyY2U6IFNpZ25hbDtcblx0X3ZhbHVlOiB1bmtub3duO1xuXHRfdmVyc2lvbjogbnVtYmVyO1xuXHRfbmV4dD86IEJhdGNoU25hcHNob3Q7XG59O1xuXG5sZXQgYmF0Y2hTbmFwc2hvdFZlcnNpb24gPSAwO1xubGV0IGN1cnJlbnRCYXRjaFNuYXBzaG90VmVyc2lvbiA9IDA7XG5sZXQgYmF0Y2hTbmFwc2hvdHM6IEJhdGNoU25hcHNob3QgfCB1bmRlZmluZWQgPSB1bmRlZmluZWQ7XG5cbi8vIEEgZ2xvYmFsIHZlcnNpb24gbnVtYmVyIGZvciBzaWduYWxzLCB1c2VkIGZvciBmYXN0LXBhdGhpbmcgcmVwZWF0ZWRcbi8vIGNvbXB1dGVkLnBlZWsoKS9jb21wdXRlZC52YWx1ZSBjYWxscyB3aGVuIG5vdGhpbmcgaGFzIGNoYW5nZWQgZ2xvYmFsbHkuXG5sZXQgZ2xvYmFsVmVyc2lvbiA9IDA7XG5cbmZ1bmN0aW9uIHJlY29yZEJhdGNoU25hcHNob3Qoc291cmNlOiBTaWduYWwpIHtcblx0Ly8gT25seSBjYXB0dXJlIHdyaXRlcyBkdXJpbmcgdGhlIHVzZXItdmlzaWJsZSBiYXRjaCBjYWxsYmFjaywgbm90IGR1cmluZyBlZmZlY3QgZmx1c2guXG5cdGlmIChiYXRjaERlcHRoID09PSAwIHx8IGJhdGNoSXRlcmF0aW9uICE9PSAwKSB7XG5cdFx0cmV0dXJuO1xuXHR9XG5cblx0aWYgKHNvdXJjZS5fYmF0Y2hTbmFwc2hvdFZlcnNpb24gIT09IGN1cnJlbnRCYXRjaFNuYXBzaG90VmVyc2lvbikge1xuXHRcdHNvdXJjZS5fYmF0Y2hTbmFwc2hvdFZlcnNpb24gPSBjdXJyZW50QmF0Y2hTbmFwc2hvdFZlcnNpb247XG5cdFx0YmF0Y2hTbmFwc2hvdHMgPSB7XG5cdFx0XHRfc291cmNlOiBzb3VyY2UsXG5cdFx0XHRfdmFsdWU6IHNvdXJjZS5fdmFsdWUsXG5cdFx0XHRfdmVyc2lvbjogc291cmNlLl92ZXJzaW9uLFxuXHRcdFx0X25leHQ6IGJhdGNoU25hcHNob3RzLFxuXHRcdH07XG5cdH1cbn1cblxuZnVuY3Rpb24gcmVjb25jaWxlQmF0Y2hTbmFwc2hvdHMoKSB7XG5cdGxldCBzbmFwc2hvdHMgPSBiYXRjaFNuYXBzaG90cztcblx0YmF0Y2hTbmFwc2hvdHMgPSB1bmRlZmluZWQ7XG5cblx0d2hpbGUgKHNuYXBzaG90cyAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0Y29uc3Qgc291cmNlID0gc25hcHNob3RzLl9zb3VyY2U7XG5cdFx0aWYgKHNvdXJjZS5fdmFsdWUgPT09IHNuYXBzaG90cy5fdmFsdWUpIHtcblx0XHRcdC8vIFRoZSB2YWx1ZSB3YXMgcmV2ZXJ0ZWQgdG8gaXRzIHByZS1iYXRjaCBzdGF0ZS4gVmVyc2lvbiBudW1iZXJzIG11c3Rcblx0XHRcdC8vIHN0YXkgbW9ub3RvbmljOiBhIGxhenkgY29tcHV0ZWQgbWF5IGhhdmUgb2JzZXJ2ZWQgYW4gaW50ZXJtZWRpYXRlXG5cdFx0XHQvLyB2ZXJzaW9uIGR1cmluZyB0aGUgYmF0Y2gsIGFuZCByb2xsaW5nIHRoZSB2ZXJzaW9uIGJhY2sgd291bGQgbGV0IGFcblx0XHRcdC8vIGZ1dHVyZSB3cml0ZSByZS1taW50IHRoYXQgb2JzZXJ2ZWQgbnVtYmVyIGZvciBhIGRpZmZlcmVudCB2YWx1ZSxcblx0XHRcdC8vIG1ha2luZyB0aGUgY29tcHV0ZWQgdHJlYXQgaXQgYXMgdW5jaGFuZ2VkIGZvcmV2ZXIuIEluc3RlYWQsXG5cdFx0XHQvLyBmYXN0LWZvcndhcmQgc3Vic2NyaWJlcnMgdGhhdCBsYXN0IHNhdyB0aGUgcHJlLWJhdGNoIHZlcnNpb24gc29cblx0XHRcdC8vIHRoZXkgc2tpcCByZWNvbXB1dGluZyBmb3IgdGhlIG5vLW9wIGNoYW5nZS5cblx0XHRcdGZvciAoXG5cdFx0XHRcdGxldCBub2RlID0gc291cmNlLl90YXJnZXRzO1xuXHRcdFx0XHRub2RlICE9PSB1bmRlZmluZWQ7XG5cdFx0XHRcdG5vZGUgPSBub2RlLl9uZXh0VGFyZ2V0XG5cdFx0XHQpIHtcblx0XHRcdFx0aWYgKG5vZGUuX3ZlcnNpb24gPT09IHNuYXBzaG90cy5fdmVyc2lvbikge1xuXHRcdFx0XHRcdG5vZGUuX3ZlcnNpb24gPSBzb3VyY2UuX3ZlcnNpb247XG5cdFx0XHRcdH1cblx0XHRcdH1cblx0XHR9XG5cdFx0c25hcHNob3RzID0gc25hcHNob3RzLl9uZXh0O1xuXHR9XG59XG5cbmZ1bmN0aW9uIGFkZERlcGVuZGVuY3koc2lnbmFsOiBTaWduYWwpOiBOb2RlIHwgdW5kZWZpbmVkIHtcblx0aWYgKGV2YWxDb250ZXh0ID09PSB1bmRlZmluZWQpIHtcblx0XHRyZXR1cm4gdW5kZWZpbmVkO1xuXHR9XG5cblx0bGV0IG5vZGUgPSBzaWduYWwuX25vZGU7XG5cdGlmIChub2RlID09PSB1bmRlZmluZWQgfHwgbm9kZS5fdGFyZ2V0ICE9PSBldmFsQ29udGV4dCkge1xuXHRcdC8qKlxuXHRcdCAqIGBzaWduYWxgIGlzIGEgbmV3IGRlcGVuZGVuY3kuIENyZWF0ZSBhIG5ldyBkZXBlbmRlbmN5IG5vZGUsIGFuZCBzZXQgaXRcblx0XHQgKiBhcyB0aGUgdGFpbCBvZiB0aGUgY3VycmVudCBjb250ZXh0J3MgZGVwZW5kZW5jeSBsaXN0LiBlLmc6XG5cdFx0ICpcblx0XHQgKiB7IEEgPC0+IEIgICAgICAgfVxuXHRcdCAqICAgICAgICAg4oaRICAgICDihpFcblx0XHQgKiAgICAgICAgdGFpbCAgbm9kZSAobmV3KVxuXHRcdCAqICAgICAgICAgICAgICAg4oaTXG5cdFx0ICogeyBBIDwtPiBCIDwtPiBDIH1cblx0XHQgKiAgICAgICAgICAgICAgIOKGkVxuXHRcdCAqICAgICAgICAgICAgICB0YWlsIChldmFsQ29udGV4dC5fc291cmNlcylcblx0XHQgKi9cblx0XHRub2RlID0ge1xuXHRcdFx0X3ZlcnNpb246IDAsXG5cdFx0XHRfc291cmNlOiBzaWduYWwsXG5cdFx0XHRfcHJldlNvdXJjZTogZXZhbENvbnRleHQuX3NvdXJjZXMsXG5cdFx0XHRfbmV4dFNvdXJjZTogdW5kZWZpbmVkLFxuXHRcdFx0X3RhcmdldDogZXZhbENvbnRleHQsXG5cdFx0XHRfcHJldlRhcmdldDogdW5kZWZpbmVkLFxuXHRcdFx0X25leHRUYXJnZXQ6IHVuZGVmaW5lZCxcblx0XHRcdF9yb2xsYmFja05vZGU6IG5vZGUsXG5cdFx0fTtcblxuXHRcdGlmIChldmFsQ29udGV4dC5fc291cmNlcyAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRldmFsQ29udGV4dC5fc291cmNlcy5fbmV4dFNvdXJjZSA9IG5vZGU7XG5cdFx0fVxuXHRcdGV2YWxDb250ZXh0Ll9zb3VyY2VzID0gbm9kZTtcblx0XHRzaWduYWwuX25vZGUgPSBub2RlO1xuXG5cdFx0Ly8gU3Vic2NyaWJlIHRvIGNoYW5nZSBub3RpZmljYXRpb25zIGZyb20gdGhpcyBkZXBlbmRlbmN5IGlmIHdlJ3JlIGluIGFuIGVmZmVjdFxuXHRcdC8vIE9SIGV2YWx1YXRpbmcgYSBjb21wdXRlZCBzaWduYWwgdGhhdCBpbiB0dXJuIGhhcyBzdWJzY3JpYmVycy5cblx0XHRpZiAoZXZhbENvbnRleHQuX2ZsYWdzICYgVFJBQ0tJTkcpIHtcblx0XHRcdHNpZ25hbC5fc3Vic2NyaWJlKG5vZGUpO1xuXHRcdH1cblx0XHRyZXR1cm4gbm9kZTtcblx0fSBlbHNlIGlmIChub2RlLl92ZXJzaW9uID09PSAtMSkge1xuXHRcdC8vIGBzaWduYWxgIGlzIGFuIGV4aXN0aW5nIGRlcGVuZGVuY3kgZnJvbSBhIHByZXZpb3VzIGV2YWx1YXRpb24uIFJldXNlIGl0LlxuXHRcdG5vZGUuX3ZlcnNpb24gPSAwO1xuXG5cdFx0LyoqXG5cdFx0ICogSWYgYG5vZGVgIGlzIG5vdCBhbHJlYWR5IHRoZSBjdXJyZW50IHRhaWwgb2YgdGhlIGRlcGVuZGVuY3kgbGlzdCAoaS5lLlxuXHRcdCAqIHRoZXJlIGlzIGEgbmV4dCBub2RlIGluIHRoZSBsaXN0KSwgdGhlbiBtYWtlIHRoZSBgbm9kZWAgdGhlIG5ldyB0YWlsLiBlLmc6XG5cdFx0ICpcblx0XHQgKiB7IEEgPC0+IEIgPC0+IEMgPC0+IEQgfVxuXHRcdCAqICAgICAgICAg4oaRICAgICAgICAgICDihpFcblx0XHQgKiAgICAgICAgbm9kZSAgIOKUjOKUgOKUgOKUgCB0YWlsIChldmFsQ29udGV4dC5fc291cmNlcylcblx0XHQgKiAgICAgICAgIOKUlOKUgOKUgOKUgOKUgOKUgOKUguKUgOKUgOKUgOKUgOKUgOKUkFxuXHRcdCAqICAgICAgICAgICAgICAg4oaTICAgICDihpNcblx0XHQgKiB7IEEgPC0+IEMgPC0+IEQgPC0+IEIgfVxuXHRcdCAqICAgICAgICAgICAgICAgICAgICAg4oaRXG5cdFx0ICogICAgICAgICAgICAgICAgICAgIHRhaWwgKGV2YWxDb250ZXh0Ll9zb3VyY2VzKVxuXHRcdCAqL1xuXHRcdGlmIChub2RlLl9uZXh0U291cmNlICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdG5vZGUuX25leHRTb3VyY2UuX3ByZXZTb3VyY2UgPSBub2RlLl9wcmV2U291cmNlO1xuXG5cdFx0XHRpZiAobm9kZS5fcHJldlNvdXJjZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRcdG5vZGUuX3ByZXZTb3VyY2UuX25leHRTb3VyY2UgPSBub2RlLl9uZXh0U291cmNlO1xuXHRcdFx0fVxuXG5cdFx0XHRub2RlLl9wcmV2U291cmNlID0gZXZhbENvbnRleHQuX3NvdXJjZXM7XG5cdFx0XHRub2RlLl9uZXh0U291cmNlID0gdW5kZWZpbmVkO1xuXG5cdFx0XHRldmFsQ29udGV4dC5fc291cmNlcyEuX25leHRTb3VyY2UgPSBub2RlO1xuXHRcdFx0ZXZhbENvbnRleHQuX3NvdXJjZXMgPSBub2RlO1xuXHRcdH1cblxuXHRcdC8vIFdlIGNhbiBhc3N1bWUgdGhhdCB0aGUgY3VycmVudGx5IGV2YWx1YXRlZCBlZmZlY3QgLyBjb21wdXRlZCBzaWduYWwgaXMgYWxyZWFkeVxuXHRcdC8vIHN1YnNjcmliZWQgdG8gY2hhbmdlIG5vdGlmaWNhdGlvbnMgZnJvbSBgc2lnbmFsYCBpZiBuZWVkZWQuXG5cdFx0cmV0dXJuIG5vZGU7XG5cdH1cblx0cmV0dXJuIHVuZGVmaW5lZDtcbn1cblxuLy8jcmVnaW9uIFNpZ25hbFxuXG4vKipcbiAqIFRoZSBiYXNlIGNsYXNzIGZvciBwbGFpbiBhbmQgY29tcHV0ZWQgc2lnbmFscy5cbiAqL1xuLy9cbi8vIEEgZnVuY3Rpb24gd2l0aCB0aGUgc2FtZSBuYW1lIGlzIGRlZmluZWQgbGF0ZXIsIHNvIHdlIG5lZWQgdG8gaWdub3JlIFR5cGVTY3JpcHQnc1xuLy8gd2FybmluZyBhYm91dCBhIHJlZGVjbGFyZWQgdmFyaWFibGUuXG4vL1xuLy8gVGhlIGNsYXNzIGlzIGRlY2xhcmVkIGhlcmUsIGJ1dCBsYXRlciBpbXBsZW1lbnRlZCB3aXRoIEVTNS1zdHlsZSBwcm90b3R5cGVzLlxuLy8gVGhpcyBlbmFibGVzIGJldHRlciBjb250cm9sIG9mIHRoZSB0cmFuc3BpbGVkIG91dHB1dCBzaXplLlxuLy8gQHRzLWlnbm9yZTogXCJDYW5ub3QgcmVkZWNsYXJlIGV4cG9ydGVkIHZhcmlhYmxlICdTaWduYWwnLlwiXG5kZWNsYXJlIGNsYXNzIFNpZ25hbDxUID0gYW55PiB7XG5cdC8qKiBAaW50ZXJuYWwgKi9cblx0X3ZhbHVlOiB1bmtub3duO1xuXG5cdC8qKlxuXHQgKiBAaW50ZXJuYWxcblx0ICogVmVyc2lvbiBudW1iZXJzIHNob3VsZCBhbHdheXMgYmUgPj0gMCwgYmVjYXVzZSB0aGUgc3BlY2lhbCB2YWx1ZSAtMSBpcyB1c2VkXG5cdCAqIGJ5IE5vZGVzIHRvIHNpZ25pZnkgcG90ZW50aWFsbHkgdW51c2VkIGJ1dCByZWN5Y2xhYmxlIG5vZGVzLlxuXHQgKi9cblx0X3ZlcnNpb246IG51bWJlcjtcblxuXHQvKiogQGludGVybmFsICovXG5cdF9ub2RlPzogTm9kZTtcblxuXHQvKiogQGludGVybmFsICovXG5cdF90YXJnZXRzPzogTm9kZTtcblxuXHQvKiogQGludGVybmFsICovXG5cdF9iYXRjaFNuYXBzaG90VmVyc2lvbjogbnVtYmVyO1xuXG5cdGNvbnN0cnVjdG9yKHZhbHVlPzogVCwgb3B0aW9ucz86IFNpZ25hbE9wdGlvbnM8VD4pO1xuXG5cdC8qKiBAaW50ZXJuYWwgKi9cblx0X3JlZnJlc2goKTogYm9vbGVhbjtcblxuXHQvKiogQGludGVybmFsICovXG5cdF9zdWJzY3JpYmUobm9kZTogTm9kZSk6IHZvaWQ7XG5cblx0LyoqIEBpbnRlcm5hbCAqL1xuXHRfdW5zdWJzY3JpYmUobm9kZTogTm9kZSk6IHZvaWQ7XG5cblx0LyoqIEBpbnRlcm5hbCAqL1xuXHRfd2F0Y2hlZD8odGhpczogU2lnbmFsPFQ+KTogdm9pZDtcblxuXHQvKiogQGludGVybmFsICovXG5cdF91bndhdGNoZWQ/KHRoaXM6IFNpZ25hbDxUPik6IHZvaWQ7XG5cblx0c3Vic2NyaWJlKGZuOiAodmFsdWU6IFQpID0+IHZvaWQpOiAoKSA9PiB2b2lkO1xuXG5cdG5hbWU/OiBzdHJpbmc7XG5cblx0dmFsdWVPZigpOiBUO1xuXG5cdHRvU3RyaW5nKCk6IHN0cmluZztcblxuXHR0b0pTT04oKTogVDtcblxuXHRwZWVrKCk6IFQ7XG5cblx0YnJhbmQ6IHR5cGVvZiBCUkFORF9TWU1CT0w7XG5cblx0Z2V0IHZhbHVlKCk6IFQ7XG5cdHNldCB2YWx1ZSh2YWx1ZTogVCk7XG59XG5cbmV4cG9ydCBpbnRlcmZhY2UgU2lnbmFsT3B0aW9uczxUID0gYW55PiB7XG5cdHdhdGNoZWQ/OiAodGhpczogU2lnbmFsPFQ+KSA9PiB2b2lkO1xuXHR1bndhdGNoZWQ/OiAodGhpczogU2lnbmFsPFQ+KSA9PiB2b2lkO1xuXHRuYW1lPzogc3RyaW5nO1xufVxuXG4vKiogQGludGVybmFsICovXG4vLyBBIGNsYXNzIHdpdGggdGhlIHNhbWUgbmFtZSBoYXMgYWxyZWFkeSBiZWVuIGRlY2xhcmVkLCBzbyB3ZSBuZWVkIHRvIGlnbm9yZVxuLy8gVHlwZVNjcmlwdCdzIHdhcm5pbmcgYWJvdXQgYSByZWRlY2xhcmVkIHZhcmlhYmxlLlxuLy9cbi8vIFRoZSBwcmV2aW91c2x5IGRlY2xhcmVkIGNsYXNzIGlzIGltcGxlbWVudGVkIGhlcmUgd2l0aCBFUzUtc3R5bGUgcHJvdG90eXBlcy5cbi8vIFRoaXMgZW5hYmxlcyBiZXR0ZXIgY29udHJvbCBvZiB0aGUgdHJhbnNwaWxlZCBvdXRwdXQgc2l6ZS5cbi8vIEB0cy1pZ25vcmU6IFwiQ2Fubm90IHJlZGVjbGFyZSBleHBvcnRlZCB2YXJpYWJsZSAnU2lnbmFsJy5cIlxuZnVuY3Rpb24gU2lnbmFsKHRoaXM6IFNpZ25hbCwgdmFsdWU/OiB1bmtub3duLCBvcHRpb25zPzogU2lnbmFsT3B0aW9ucykge1xuXHR0aGlzLl92YWx1ZSA9IHZhbHVlO1xuXHR0aGlzLl92ZXJzaW9uID0gMDtcblx0dGhpcy5fbm9kZSA9IHVuZGVmaW5lZDtcblx0dGhpcy5fdGFyZ2V0cyA9IHVuZGVmaW5lZDtcblx0dGhpcy5fYmF0Y2hTbmFwc2hvdFZlcnNpb24gPSAwO1xuXHR0aGlzLl93YXRjaGVkID0gb3B0aW9ucz8ud2F0Y2hlZDtcblx0dGhpcy5fdW53YXRjaGVkID0gb3B0aW9ucz8udW53YXRjaGVkO1xuXHR0aGlzLm5hbWUgPSBvcHRpb25zPy5uYW1lO1xufVxuXG5TaWduYWwucHJvdG90eXBlLmJyYW5kID0gQlJBTkRfU1lNQk9MO1xuXG5TaWduYWwucHJvdG90eXBlLl9yZWZyZXNoID0gZnVuY3Rpb24gKCkge1xuXHRyZXR1cm4gdHJ1ZTtcbn07XG5cblNpZ25hbC5wcm90b3R5cGUuX3N1YnNjcmliZSA9IGZ1bmN0aW9uIChub2RlKSB7XG5cdGNvbnN0IHRhcmdldHMgPSB0aGlzLl90YXJnZXRzO1xuXHRpZiAodGFyZ2V0cyAhPT0gbm9kZSAmJiBub2RlLl9wcmV2VGFyZ2V0ID09PSB1bmRlZmluZWQpIHtcblx0XHRub2RlLl9uZXh0VGFyZ2V0ID0gdGFyZ2V0cztcblx0XHR0aGlzLl90YXJnZXRzID0gbm9kZTtcblxuXHRcdGlmICh0YXJnZXRzICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdHRhcmdldHMuX3ByZXZUYXJnZXQgPSBub2RlO1xuXHRcdH0gZWxzZSB7XG5cdFx0XHR1bnRyYWNrZWQoKCkgPT4ge1xuXHRcdFx0XHR0aGlzLl93YXRjaGVkPy5jYWxsKHRoaXMpO1xuXHRcdFx0fSk7XG5cdFx0fVxuXHR9XG59O1xuXG5TaWduYWwucHJvdG90eXBlLl91bnN1YnNjcmliZSA9IGZ1bmN0aW9uIChub2RlKSB7XG5cdC8vIE9ubHkgcnVuIHRoZSB1bnN1YnNjcmliZSBzdGVwIGlmIHRoZSBzaWduYWwgaGFzIGFueSBzdWJzY3JpYmVycyB0byBiZWdpbiB3aXRoLlxuXHRpZiAodGhpcy5fdGFyZ2V0cyAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0Y29uc3QgcHJldiA9IG5vZGUuX3ByZXZUYXJnZXQ7XG5cdFx0Y29uc3QgbmV4dCA9IG5vZGUuX25leHRUYXJnZXQ7XG5cdFx0aWYgKHByZXYgIT09IHVuZGVmaW5lZCkge1xuXHRcdFx0cHJldi5fbmV4dFRhcmdldCA9IG5leHQ7XG5cdFx0XHRub2RlLl9wcmV2VGFyZ2V0ID0gdW5kZWZpbmVkO1xuXHRcdH1cblxuXHRcdGlmIChuZXh0ICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdG5leHQuX3ByZXZUYXJnZXQgPSBwcmV2O1xuXHRcdFx0bm9kZS5fbmV4dFRhcmdldCA9IHVuZGVmaW5lZDtcblx0XHR9XG5cblx0XHRpZiAobm9kZSA9PT0gdGhpcy5fdGFyZ2V0cykge1xuXHRcdFx0dGhpcy5fdGFyZ2V0cyA9IG5leHQ7XG5cdFx0XHRpZiAobmV4dCA9PT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRcdHVudHJhY2tlZCgoKSA9PiB7XG5cdFx0XHRcdFx0dGhpcy5fdW53YXRjaGVkPy5jYWxsKHRoaXMpO1xuXHRcdFx0XHR9KTtcblx0XHRcdH1cblx0XHR9XG5cdH1cbn07XG5cblNpZ25hbC5wcm90b3R5cGUuc3Vic2NyaWJlID0gZnVuY3Rpb24gKGZuKSB7XG5cdHJldHVybiBlZmZlY3QoXG5cdFx0KCkgPT4ge1xuXHRcdFx0Y29uc3QgdmFsdWUgPSB0aGlzLnZhbHVlO1xuXHRcdFx0dW50cmFja2VkKCgpID0+IGZuKHZhbHVlKSk7XG5cdFx0fSxcblx0XHR7IG5hbWU6IFwic3ViXCIgfVxuXHQpO1xufTtcblxuU2lnbmFsLnByb3RvdHlwZS52YWx1ZU9mID0gZnVuY3Rpb24gKCkge1xuXHRyZXR1cm4gdGhpcy52YWx1ZTtcbn07XG5cblNpZ25hbC5wcm90b3R5cGUudG9TdHJpbmcgPSBmdW5jdGlvbiAoKSB7XG5cdHJldHVybiB0aGlzLnZhbHVlICsgXCJcIjtcbn07XG5cblNpZ25hbC5wcm90b3R5cGUudG9KU09OID0gZnVuY3Rpb24gKCkge1xuXHRyZXR1cm4gdGhpcy52YWx1ZTtcbn07XG5cblNpZ25hbC5wcm90b3R5cGUucGVlayA9IGZ1bmN0aW9uICgpIHtcblx0cmV0dXJuIHVudHJhY2tlZCgoKSA9PiB0aGlzLnZhbHVlKTtcbn07XG5cbk9iamVjdC5kZWZpbmVQcm9wZXJ0eShTaWduYWwucHJvdG90eXBlLCBcInZhbHVlXCIsIHtcblx0Z2V0KHRoaXM6IFNpZ25hbCkge1xuXHRcdGNvbnN0IG5vZGUgPSBhZGREZXBlbmRlbmN5KHRoaXMpO1xuXHRcdGlmIChub2RlICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdG5vZGUuX3ZlcnNpb24gPSB0aGlzLl92ZXJzaW9uO1xuXHRcdH1cblx0XHRyZXR1cm4gdGhpcy5fdmFsdWU7XG5cdH0sXG5cdHNldCh0aGlzOiBTaWduYWwsIHZhbHVlKSB7XG5cdFx0aWYgKHZhbHVlICE9PSB0aGlzLl92YWx1ZSkge1xuXHRcdFx0aWYgKGJhdGNoSXRlcmF0aW9uID4gMTAwKSB7XG5cdFx0XHRcdHRocm93IG5ldyBFcnJvcihcIkN5Y2xlIGRldGVjdGVkXCIpO1xuXHRcdFx0fVxuXG5cdFx0XHRyZWNvcmRCYXRjaFNuYXBzaG90KHRoaXMpO1xuXHRcdFx0dGhpcy5fdmFsdWUgPSB2YWx1ZTtcblx0XHRcdHRoaXMuX3ZlcnNpb24rKztcblx0XHRcdGdsb2JhbFZlcnNpb24rKztcblxuXHRcdFx0LyoqQF9fSU5MSU5FX18qLyBzdGFydEJhdGNoKCk7XG5cdFx0XHR0cnkge1xuXHRcdFx0XHRmb3IgKFxuXHRcdFx0XHRcdGxldCBub2RlID0gdGhpcy5fdGFyZ2V0cztcblx0XHRcdFx0XHRub2RlICE9PSB1bmRlZmluZWQ7XG5cdFx0XHRcdFx0bm9kZSA9IG5vZGUuX25leHRUYXJnZXRcblx0XHRcdFx0KSB7XG5cdFx0XHRcdFx0bm9kZS5fdGFyZ2V0Ll9ub3RpZnkoKTtcblx0XHRcdFx0fVxuXHRcdFx0fSBmaW5hbGx5IHtcblx0XHRcdFx0ZW5kQmF0Y2goKTtcblx0XHRcdH1cblx0XHR9XG5cdH0sXG59KTtcblxuLyoqXG4gKiBDcmVhdGUgYSBuZXcgcGxhaW4gc2lnbmFsLlxuICpcbiAqIEBwYXJhbSB2YWx1ZSBUaGUgaW5pdGlhbCB2YWx1ZSBmb3IgdGhlIHNpZ25hbC5cbiAqIEByZXR1cm5zIEEgbmV3IHNpZ25hbC5cbiAqL1xuZXhwb3J0IGZ1bmN0aW9uIHNpZ25hbDxUPih2YWx1ZTogVCwgb3B0aW9ucz86IFNpZ25hbE9wdGlvbnM8VD4pOiBTaWduYWw8VD47XG5leHBvcnQgZnVuY3Rpb24gc2lnbmFsPFQgPSB1bmRlZmluZWQ+KCk6IFNpZ25hbDxUIHwgdW5kZWZpbmVkPjtcbmV4cG9ydCBmdW5jdGlvbiBzaWduYWw8VD4odmFsdWU/OiBULCBvcHRpb25zPzogU2lnbmFsT3B0aW9uczxUPik6IFNpZ25hbDxUPiB7XG5cdHJldHVybiBuZXcgU2lnbmFsKHZhbHVlLCBvcHRpb25zKTtcbn1cblxuLy8jZW5kcmVnaW9uIFNpZ25hbFxuXG4vLyNyZWdpb24gQ29tcHV0ZWRcblxuZnVuY3Rpb24gbmVlZHNUb1JlY29tcHV0ZSh0YXJnZXQ6IENvbXB1dGVkIHwgRWZmZWN0KTogYm9vbGVhbiB7XG5cdC8vIENoZWNrIHRoZSBkZXBlbmRlbmNpZXMgZm9yIGNoYW5nZWQgdmFsdWVzLiBUaGUgZGVwZW5kZW5jeSBsaXN0IGlzIGFscmVhZHlcblx0Ly8gaW4gb3JkZXIgb2YgdXNlLiBUaGVyZWZvcmUgaWYgbXVsdGlwbGUgZGVwZW5kZW5jaWVzIGhhdmUgY2hhbmdlZCB2YWx1ZXMsIG9ubHlcblx0Ly8gdGhlIGZpcnN0IHVzZWQgZGVwZW5kZW5jeSBpcyByZS1ldmFsdWF0ZWQgYXQgdGhpcyBwb2ludC5cblx0Zm9yIChcblx0XHRsZXQgbm9kZSA9IHRhcmdldC5fc291cmNlcztcblx0XHRub2RlICE9PSB1bmRlZmluZWQ7XG5cdFx0bm9kZSA9IG5vZGUuX25leHRTb3VyY2Vcblx0KSB7XG5cdFx0aWYgKFxuXHRcdFx0Ly8gSWYgdGhlIGRlcGVuZGVuY3kgaGFzIGRlZmluaXRlbHkgYmVlbiB1cGRhdGVkIHNpbmNlIGl0cyB2ZXJzaW9uIG51bWJlclxuXHRcdFx0Ly8gd2FzIG9ic2VydmVkLCB0aGVuIHdlIG5lZWQgdG8gcmVjb21wdXRlLiBUaGlzIGZpcnN0IGNoZWNrIGlzIG5vdCBzdHJpY3RseVxuXHRcdFx0Ly8gbmVjZXNzYXJ5IGZvciBjb3JyZWN0bmVzcywgYnV0IGFsbG93cyB1cyB0byBza2lwIHRoZSByZWZyZXNoIGNhbGwgaWYgdGhlXG5cdFx0XHQvLyBkZXBlbmRlbmN5IGhhcyBhbHJlYWR5IGJlZW4gdXBkYXRlZC5cblx0XHRcdG5vZGUuX3NvdXJjZS5fdmVyc2lvbiAhPT0gbm9kZS5fdmVyc2lvbiB8fFxuXHRcdFx0Ly8gUmVmcmVzaCB0aGUgZGVwZW5kZW5jeS4gSWYgdGhlcmUncyBzb21ldGhpbmcgYmxvY2tpbmcgdGhlIHJlZnJlc2ggKGUuZy4gYVxuXHRcdFx0Ly8gZGVwZW5kZW5jeSBjeWNsZSksIHRoZW4gd2UgbmVlZCB0byByZWNvbXB1dGUuXG5cdFx0XHQhbm9kZS5fc291cmNlLl9yZWZyZXNoKCkgfHxcblx0XHRcdC8vIElmIHRoZSBkZXBlbmRlbmN5IGdvdCBhIG5ldyB2ZXJzaW9uIGFmdGVyIHRoZSByZWZyZXNoLCB0aGVuIHdlIG5lZWQgdG8gcmVjb21wdXRlLlxuXHRcdFx0bm9kZS5fc291cmNlLl92ZXJzaW9uICE9PSBub2RlLl92ZXJzaW9uXG5cdFx0KSB7XG5cdFx0XHRyZXR1cm4gdHJ1ZTtcblx0XHR9XG5cdH1cblx0Ly8gSWYgbm9uZSBvZiB0aGUgZGVwZW5kZW5jaWVzIGhhdmUgY2hhbmdlZCB2YWx1ZXMgc2luY2UgbGFzdCByZWNvbXB1dGUgdGhlblxuXHQvLyB0aGVyZSdzIG5vIG5lZWQgdG8gcmVjb21wdXRlLlxuXHRyZXR1cm4gZmFsc2U7XG59XG5cbmZ1bmN0aW9uIHByZXBhcmVTb3VyY2VzKHRhcmdldDogQ29tcHV0ZWQgfCBFZmZlY3QpIHtcblx0LyoqXG5cdCAqIDEuIE1hcmsgYWxsIGN1cnJlbnQgc291cmNlcyBhcyByZS11c2FibGUgbm9kZXMgKHZlcnNpb246IC0xKVxuXHQgKiAyLiBTZXQgYSByb2xsYmFjayBub2RlIGlmIHRoZSBjdXJyZW50IG5vZGUgaXMgYmVpbmcgdXNlZCBpbiBhIGRpZmZlcmVudCBjb250ZXh0XG5cdCAqIDMuIFBvaW50ICd0YXJnZXQuX3NvdXJjZXMnIHRvIHRoZSB0YWlsIG9mIHRoZSBkb3VibHktbGlua2VkIGxpc3QsIGUuZzpcblx0ICpcblx0ICogICAgeyB1bmRlZmluZWQgPC0gQSA8LT4gQiA8LT4gQyAtPiB1bmRlZmluZWQgfVxuXHQgKiAgICAgICAgICAgICAgICAgICDihpEgICAgICAgICAgIOKGkVxuXHQgKiAgICAgICAgICAgICAgICAgICDilIIgICAgICAgICAgIOKUlOKUgOKUgOKUgOKUgOKUgOKUgOKUkFxuXHQgKiB0YXJnZXQuX3NvdXJjZXMgPSBBOyAobm9kZSBpcyBoZWFkKSAg4pSCXG5cdCAqICAgICAgICAgICAgICAgICAgIOKGkyAgICAgICAgICAgICAgICAgIOKUglxuXHQgKiB0YXJnZXQuX3NvdXJjZXMgPSBDOyAobm9kZSBpcyB0YWlsKSDilIDilJhcblx0ICovXG5cdGZvciAoXG5cdFx0bGV0IG5vZGUgPSB0YXJnZXQuX3NvdXJjZXM7XG5cdFx0bm9kZSAhPT0gdW5kZWZpbmVkO1xuXHRcdG5vZGUgPSBub2RlLl9uZXh0U291cmNlXG5cdCkge1xuXHRcdGNvbnN0IHJvbGxiYWNrTm9kZSA9IG5vZGUuX3NvdXJjZS5fbm9kZTtcblx0XHRpZiAocm9sbGJhY2tOb2RlICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdG5vZGUuX3JvbGxiYWNrTm9kZSA9IHJvbGxiYWNrTm9kZTtcblx0XHR9XG5cdFx0bm9kZS5fc291cmNlLl9ub2RlID0gbm9kZTtcblx0XHRub2RlLl92ZXJzaW9uID0gLTE7XG5cblx0XHRpZiAobm9kZS5fbmV4dFNvdXJjZSA9PT0gdW5kZWZpbmVkKSB7XG5cdFx0XHR0YXJnZXQuX3NvdXJjZXMgPSBub2RlO1xuXHRcdFx0YnJlYWs7XG5cdFx0fVxuXHR9XG59XG5cbmZ1bmN0aW9uIGNsZWFudXBTb3VyY2VzKHRhcmdldDogQ29tcHV0ZWQgfCBFZmZlY3QpIHtcblx0bGV0IG5vZGUgPSB0YXJnZXQuX3NvdXJjZXM7XG5cdGxldCBoZWFkOiBOb2RlIHwgdW5kZWZpbmVkID0gdW5kZWZpbmVkO1xuXG5cdC8qKlxuXHQgKiBBdCB0aGlzIHBvaW50ICd0YXJnZXQuX3NvdXJjZXMnIHBvaW50cyB0byB0aGUgdGFpbCBvZiB0aGUgZG91Ymx5LWxpbmtlZCBsaXN0LlxuXHQgKiBJdCBjb250YWlucyBhbGwgZXhpc3Rpbmcgc291cmNlcyArIG5ldyBzb3VyY2VzIGluIG9yZGVyIG9mIHVzZS5cblx0ICogSXRlcmF0ZSBiYWNrd2FyZHMgdW50aWwgd2UgZmluZCB0aGUgaGVhZCBub2RlIHdoaWxlIGRyb3BwaW5nIG9sZCBkZXBlbmRlbmNpZXMuXG5cdCAqL1xuXHR3aGlsZSAobm9kZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0Y29uc3QgcHJldiA9IG5vZGUuX3ByZXZTb3VyY2U7XG5cblx0XHQvKipcblx0XHQgKiBUaGUgbm9kZSB3YXMgbm90IHJlLXVzZWQsIHVuc3Vic2NyaWJlIGZyb20gaXRzIGNoYW5nZSBub3RpZmljYXRpb25zIGFuZCByZW1vdmUgaXRzZWxmXG5cdFx0ICogZnJvbSB0aGUgZG91Ymx5LWxpbmtlZCBsaXN0LiBlLmc6XG5cdFx0ICpcblx0XHQgKiB7IEEgPC0+IEIgPC0+IEMgfVxuXHRcdCAqICAgICAgICAg4oaTXG5cdFx0ICogICAgeyBBIDwtPiBDIH1cblx0XHQgKi9cblx0XHRpZiAobm9kZS5fdmVyc2lvbiA9PT0gLTEpIHtcblx0XHRcdG5vZGUuX3NvdXJjZS5fdW5zdWJzY3JpYmUobm9kZSk7XG5cblx0XHRcdGlmIChwcmV2ICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdFx0cHJldi5fbmV4dFNvdXJjZSA9IG5vZGUuX25leHRTb3VyY2U7XG5cdFx0XHR9XG5cdFx0XHRpZiAobm9kZS5fbmV4dFNvdXJjZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRcdG5vZGUuX25leHRTb3VyY2UuX3ByZXZTb3VyY2UgPSBwcmV2O1xuXHRcdFx0fVxuXHRcdH0gZWxzZSB7XG5cdFx0XHQvKipcblx0XHRcdCAqIFRoZSBuZXcgaGVhZCBpcyB0aGUgbGFzdCBub2RlIHNlZW4gd2hpY2ggd2Fzbid0IHJlbW92ZWQvdW5zdWJzY3JpYmVkXG5cdFx0XHQgKiBmcm9tIHRoZSBkb3VibHktbGlua2VkIGxpc3QuIGUuZzpcblx0XHRcdCAqXG5cdFx0XHQgKiB7IEEgPC0+IEIgPC0+IEMgfVxuXHRcdFx0ICogICDihpEgICAgIOKGkSAgICAg4oaRXG5cdFx0XHQgKiAgIOKUgiAgICAg4pSCICAgICDilJQgaGVhZCA9IG5vZGVcblx0XHRcdCAqICAg4pSCICAgICDilJQgaGVhZCA9IG5vZGVcblx0XHRcdCAqICAg4pSUIGhlYWQgPSBub2RlXG5cdFx0XHQgKi9cblx0XHRcdGhlYWQgPSBub2RlO1xuXHRcdH1cblxuXHRcdG5vZGUuX3NvdXJjZS5fbm9kZSA9IG5vZGUuX3JvbGxiYWNrTm9kZTtcblx0XHRpZiAobm9kZS5fcm9sbGJhY2tOb2RlICE9PSB1bmRlZmluZWQpIHtcblx0XHRcdG5vZGUuX3JvbGxiYWNrTm9kZSA9IHVuZGVmaW5lZDtcblx0XHR9XG5cblx0XHRub2RlID0gcHJldjtcblx0fVxuXG5cdHRhcmdldC5fc291cmNlcyA9IGhlYWQ7XG59XG5cbi8qKlxuICogVGhlIGJhc2UgY2xhc3MgZm9yIGNvbXB1dGVkIHNpZ25hbHMuXG4gKi9cbmRlY2xhcmUgY2xhc3MgQ29tcHV0ZWQ8VCA9IGFueT4gZXh0ZW5kcyBTaWduYWw8VD4ge1xuXHRfZm46ICgpID0+IFQ7XG5cdF9zb3VyY2VzPzogTm9kZTtcblx0X2dsb2JhbFZlcnNpb246IG51bWJlcjtcblx0X2ZsYWdzOiBudW1iZXI7XG5cblx0Y29uc3RydWN0b3IoZm46ICgpID0+IFQsIG9wdGlvbnM/OiBTaWduYWxPcHRpb25zPFQ+KTtcblxuXHRfbm90aWZ5KCk6IHZvaWQ7XG5cdGdldCB2YWx1ZSgpOiBUO1xufVxuXG4vKiogQGludGVybmFsICovXG5mdW5jdGlvbiBDb21wdXRlZCh0aGlzOiBDb21wdXRlZCwgZm46ICgpID0+IHVua25vd24sIG9wdGlvbnM/OiBTaWduYWxPcHRpb25zKSB7XG5cdFNpZ25hbC5jYWxsKHRoaXMsIHVuZGVmaW5lZCwgb3B0aW9ucyk7XG5cblx0dGhpcy5fZm4gPSBmbjtcblx0dGhpcy5fc291cmNlcyA9IHVuZGVmaW5lZDtcblx0dGhpcy5fZ2xvYmFsVmVyc2lvbiA9IGdsb2JhbFZlcnNpb24gLSAxO1xuXHR0aGlzLl9mbGFncyA9IE9VVERBVEVEO1xufVxuXG5Db21wdXRlZC5wcm90b3R5cGUgPSBuZXcgU2lnbmFsKCkgYXMgQ29tcHV0ZWQ7XG5cbkNvbXB1dGVkLnByb3RvdHlwZS5fcmVmcmVzaCA9IGZ1bmN0aW9uICgpIHtcblx0dGhpcy5fZmxhZ3MgJj0gfk5PVElGSUVEO1xuXG5cdGlmICh0aGlzLl9mbGFncyAmIFJVTk5JTkcpIHtcblx0XHRyZXR1cm4gZmFsc2U7XG5cdH1cblxuXHQvLyBJZiB0aGlzIGNvbXB1dGVkIHNpZ25hbCBoYXMgc3Vic2NyaWJlZCB0byB1cGRhdGVzIGZyb20gaXRzIGRlcGVuZGVuY2llc1xuXHQvLyAoVFJBQ0tJTkcgZmxhZyBzZXQpIGFuZCBub25lIG9mIHRoZW0gaGF2ZSBub3RpZmllZCBhYm91dCBjaGFuZ2VzIChPVVREQVRFRFxuXHQvLyBmbGFnIG5vdCBzZXQpLCB0aGVuIHRoZSBjb21wdXRlZCB2YWx1ZSBjYW4ndCBoYXZlIGNoYW5nZWQuXG5cdGlmICgodGhpcy5fZmxhZ3MgJiAoT1VUREFURUQgfCBUUkFDS0lORykpID09PSBUUkFDS0lORykge1xuXHRcdHJldHVybiB0cnVlO1xuXHR9XG5cdHRoaXMuX2ZsYWdzICY9IH5PVVREQVRFRDtcblxuXHRpZiAodGhpcy5fZ2xvYmFsVmVyc2lvbiA9PT0gZ2xvYmFsVmVyc2lvbikge1xuXHRcdHJldHVybiB0cnVlO1xuXHR9XG5cdHRoaXMuX2dsb2JhbFZlcnNpb24gPSBnbG9iYWxWZXJzaW9uO1xuXG5cdC8vIE1hcmsgdGhpcyBjb21wdXRlZCBzaWduYWwgcnVubmluZyBiZWZvcmUgY2hlY2tpbmcgdGhlIGRlcGVuZGVuY2llcyBmb3IgdmFsdWVcblx0Ly8gY2hhbmdlcywgc28gdGhhdCB0aGUgUlVOTklORyBmbGFnIGNhbiBiZSB1c2VkIHRvIG5vdGljZSBjeWNsaWNhbCBkZXBlbmRlbmNpZXMuXG5cdHRoaXMuX2ZsYWdzIHw9IFJVTk5JTkc7XG5cdGlmICh0aGlzLl92ZXJzaW9uID4gMCAmJiAhbmVlZHNUb1JlY29tcHV0ZSh0aGlzKSkge1xuXHRcdHRoaXMuX2ZsYWdzICY9IH5SVU5OSU5HO1xuXHRcdHJldHVybiB0cnVlO1xuXHR9XG5cblx0Y29uc3QgcHJldkNvbnRleHQgPSBldmFsQ29udGV4dDtcblx0dHJ5IHtcblx0XHRwcmVwYXJlU291cmNlcyh0aGlzKTtcblx0XHRldmFsQ29udGV4dCA9IHRoaXM7XG5cdFx0Y29uc3QgdmFsdWUgPSB0aGlzLl9mbigpO1xuXHRcdGlmIChcblx0XHRcdHRoaXMuX2ZsYWdzICYgSEFTX0VSUk9SIHx8XG5cdFx0XHR0aGlzLl92YWx1ZSAhPT0gdmFsdWUgfHxcblx0XHRcdHRoaXMuX3ZlcnNpb24gPT09IDBcblx0XHQpIHtcblx0XHRcdHRoaXMuX3ZhbHVlID0gdmFsdWU7XG5cdFx0XHR0aGlzLl9mbGFncyAmPSB+SEFTX0VSUk9SO1xuXHRcdFx0dGhpcy5fdmVyc2lvbisrO1xuXHRcdH1cblx0fSBjYXRjaCAoZXJyKSB7XG5cdFx0dGhpcy5fdmFsdWUgPSBlcnI7XG5cdFx0dGhpcy5fZmxhZ3MgfD0gSEFTX0VSUk9SO1xuXHRcdHRoaXMuX3ZlcnNpb24rKztcblx0fVxuXHRldmFsQ29udGV4dCA9IHByZXZDb250ZXh0O1xuXHRjbGVhbnVwU291cmNlcyh0aGlzKTtcblx0dGhpcy5fZmxhZ3MgJj0gflJVTk5JTkc7XG5cdHJldHVybiB0cnVlO1xufTtcblxuQ29tcHV0ZWQucHJvdG90eXBlLl9zdWJzY3JpYmUgPSBmdW5jdGlvbiAobm9kZSkge1xuXHRpZiAodGhpcy5fdGFyZ2V0cyA9PT0gdW5kZWZpbmVkKSB7XG5cdFx0dGhpcy5fZmxhZ3MgfD0gT1VUREFURUQgfCBUUkFDS0lORztcblxuXHRcdC8vIEEgY29tcHV0ZWQgc2lnbmFsIHN1YnNjcmliZXMgbGF6aWx5IHRvIGl0cyBkZXBlbmRlbmNpZXMgd2hlbiBpdFxuXHRcdC8vIGdldHMgaXRzIGZpcnN0IHN1YnNjcmliZXIuXG5cdFx0Zm9yIChcblx0XHRcdGxldCBub2RlID0gdGhpcy5fc291cmNlcztcblx0XHRcdG5vZGUgIT09IHVuZGVmaW5lZDtcblx0XHRcdG5vZGUgPSBub2RlLl9uZXh0U291cmNlXG5cdFx0KSB7XG5cdFx0XHRub2RlLl9zb3VyY2UuX3N1YnNjcmliZShub2RlKTtcblx0XHR9XG5cdH1cblx0U2lnbmFsLnByb3RvdHlwZS5fc3Vic2NyaWJlLmNhbGwodGhpcywgbm9kZSk7XG59O1xuXG5Db21wdXRlZC5wcm90b3R5cGUuX3Vuc3Vic2NyaWJlID0gZnVuY3Rpb24gKG5vZGUpIHtcblx0Ly8gT25seSBydW4gdGhlIHVuc3Vic2NyaWJlIHN0ZXAgaWYgdGhlIGNvbXB1dGVkIHNpZ25hbCBoYXMgYW55IHN1YnNjcmliZXJzLlxuXHRpZiAodGhpcy5fdGFyZ2V0cyAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0U2lnbmFsLnByb3RvdHlwZS5fdW5zdWJzY3JpYmUuY2FsbCh0aGlzLCBub2RlKTtcblxuXHRcdC8vIENvbXB1dGVkIHNpZ25hbCB1bnN1YnNjcmliZXMgZnJvbSBpdHMgZGVwZW5kZW5jaWVzIHdoZW4gaXQgbG9zZXMgaXRzIGxhc3Qgc3Vic2NyaWJlci5cblx0XHQvLyBUaGlzIG1ha2VzIGl0IHBvc3NpYmxlIGZvciB1bnJlZmVyZW5jZXMgc3ViZ3JhcGhzIG9mIGNvbXB1dGVkIHNpZ25hbHMgdG8gZ2V0IGdhcmJhZ2UgY29sbGVjdGVkLlxuXHRcdGlmICh0aGlzLl90YXJnZXRzID09PSB1bmRlZmluZWQpIHtcblx0XHRcdHRoaXMuX2ZsYWdzICY9IH5UUkFDS0lORztcblxuXHRcdFx0Zm9yIChcblx0XHRcdFx0bGV0IG5vZGUgPSB0aGlzLl9zb3VyY2VzO1xuXHRcdFx0XHRub2RlICE9PSB1bmRlZmluZWQ7XG5cdFx0XHRcdG5vZGUgPSBub2RlLl9uZXh0U291cmNlXG5cdFx0XHQpIHtcblx0XHRcdFx0bm9kZS5fc291cmNlLl91bnN1YnNjcmliZShub2RlKTtcblx0XHRcdH1cblx0XHR9XG5cdH1cbn07XG5cbkNvbXB1dGVkLnByb3RvdHlwZS5fbm90aWZ5ID0gZnVuY3Rpb24gKCkge1xuXHRpZiAoISh0aGlzLl9mbGFncyAmIE5PVElGSUVEKSkge1xuXHRcdHRoaXMuX2ZsYWdzIHw9IE9VVERBVEVEIHwgTk9USUZJRUQ7XG5cblx0XHRmb3IgKFxuXHRcdFx0bGV0IG5vZGUgPSB0aGlzLl90YXJnZXRzO1xuXHRcdFx0bm9kZSAhPT0gdW5kZWZpbmVkO1xuXHRcdFx0bm9kZSA9IG5vZGUuX25leHRUYXJnZXRcblx0XHQpIHtcblx0XHRcdG5vZGUuX3RhcmdldC5fbm90aWZ5KCk7XG5cdFx0fVxuXHR9XG59O1xuXG5PYmplY3QuZGVmaW5lUHJvcGVydHkoQ29tcHV0ZWQucHJvdG90eXBlLCBcInZhbHVlXCIsIHtcblx0Z2V0KHRoaXM6IENvbXB1dGVkKSB7XG5cdFx0aWYgKHRoaXMuX2ZsYWdzICYgUlVOTklORykge1xuXHRcdFx0dGhyb3cgbmV3IEVycm9yKFwiQ3ljbGUgZGV0ZWN0ZWRcIik7XG5cdFx0fVxuXHRcdGNvbnN0IG5vZGUgPSBhZGREZXBlbmRlbmN5KHRoaXMpO1xuXHRcdHRoaXMuX3JlZnJlc2goKTtcblx0XHRpZiAobm9kZSAhPT0gdW5kZWZpbmVkKSB7XG5cdFx0XHRub2RlLl92ZXJzaW9uID0gdGhpcy5fdmVyc2lvbjtcblx0XHR9XG5cdFx0aWYgKHRoaXMuX2ZsYWdzICYgSEFTX0VSUk9SKSB7XG5cdFx0XHR0aHJvdyB0aGlzLl92YWx1ZTtcblx0XHR9XG5cdFx0cmV0dXJuIHRoaXMuX3ZhbHVlO1xuXHR9LFxufSk7XG5cbi8qKlxuICogQW4gaW50ZXJmYWNlIGZvciByZWFkLW9ubHkgc2lnbmFscy5cbiAqL1xuaW50ZXJmYWNlIFJlYWRvbmx5U2lnbmFsPFQgPSBhbnk+IHtcblx0cmVhZG9ubHkgdmFsdWU6IFQ7XG5cdHBlZWsoKTogVDtcblxuXHRzdWJzY3JpYmUoZm46ICh2YWx1ZTogVCkgPT4gdm9pZCk6ICgpID0+IHZvaWQ7XG5cdHZhbHVlT2YoKTogVDtcblx0dG9TdHJpbmcoKTogc3RyaW5nO1xuXHR0b0pTT04oKTogVDtcblx0YnJhbmQ6IHR5cGVvZiBCUkFORF9TWU1CT0w7XG59XG5cbi8qKlxuICogQ3JlYXRlIGEgbmV3IHNpZ25hbCB0aGF0IGlzIGNvbXB1dGVkIGJhc2VkIG9uIHRoZSB2YWx1ZXMgb2Ygb3RoZXIgc2lnbmFscy5cbiAqXG4gKiBUaGUgcmV0dXJuZWQgY29tcHV0ZWQgc2lnbmFsIGlzIHJlYWQtb25seSwgYW5kIGl0cyB2YWx1ZSBpcyBhdXRvbWF0aWNhbGx5XG4gKiB1cGRhdGVkIHdoZW4gYW55IHNpZ25hbHMgYWNjZXNzZWQgZnJvbSB3aXRoaW4gdGhlIGNhbGxiYWNrIGZ1bmN0aW9uIGNoYW5nZS5cbiAqXG4gKiBAcGFyYW0gZm4gVGhlIGVmZmVjdCBjYWxsYmFjay5cbiAqIEByZXR1cm5zIEEgbmV3IHJlYWQtb25seSBzaWduYWwuXG4gKi9cbmZ1bmN0aW9uIGNvbXB1dGVkPFQ+KFxuXHRmbjogKCkgPT4gVCxcblx0b3B0aW9ucz86IFNpZ25hbE9wdGlvbnM8VD5cbik6IFJlYWRvbmx5U2lnbmFsPFQ+IHtcblx0cmV0dXJuIG5ldyBDb21wdXRlZChmbiwgb3B0aW9ucyk7XG59XG5cbi8vI2VuZHJlZ2lvbiBDb21wdXRlZFxuXG4vLyNyZWdpb24gRWZmZWN0XG5cbmZ1bmN0aW9uIGNsZWFudXBFZmZlY3QoZWZmZWN0OiBFZmZlY3QpIHtcblx0Y29uc3QgY2xlYW51cCA9IGVmZmVjdC5fY2xlYW51cDtcblx0ZWZmZWN0Ll9jbGVhbnVwID0gdW5kZWZpbmVkO1xuXG5cdGlmICh0eXBlb2YgY2xlYW51cCA9PT0gXCJmdW5jdGlvblwiKSB7XG5cdFx0LypAX19JTkxJTkVfXyoqLyBzdGFydEJhdGNoKCk7XG5cblx0XHQvLyBSdW4gY2xlYW51cCBmdW5jdGlvbnMgYWx3YXlzIG91dHNpZGUgb2YgYW55IGNvbnRleHQuXG5cdFx0Y29uc3QgcHJldkNvbnRleHQgPSBldmFsQ29udGV4dDtcblx0XHRldmFsQ29udGV4dCA9IHVuZGVmaW5lZDtcblx0XHR0cnkge1xuXHRcdFx0Y2xlYW51cCgpO1xuXHRcdH0gY2F0Y2ggKGVycikge1xuXHRcdFx0ZWZmZWN0Ll9mbGFncyAmPSB+UlVOTklORztcblx0XHRcdGVmZmVjdC5fZmxhZ3MgfD0gRElTUE9TRUQ7XG5cdFx0XHRkaXNwb3NlRWZmZWN0KGVmZmVjdCk7XG5cdFx0XHR0aHJvdyBlcnI7XG5cdFx0fSBmaW5hbGx5IHtcblx0XHRcdGV2YWxDb250ZXh0ID0gcHJldkNvbnRleHQ7XG5cdFx0XHRlbmRCYXRjaCgpO1xuXHRcdH1cblx0fVxufVxuXG5mdW5jdGlvbiBkaXNwb3NlRWZmZWN0KGVmZmVjdDogRWZmZWN0KSB7XG5cdGZvciAoXG5cdFx0bGV0IG5vZGUgPSBlZmZlY3QuX3NvdXJjZXM7XG5cdFx0bm9kZSAhPT0gdW5kZWZpbmVkO1xuXHRcdG5vZGUgPSBub2RlLl9uZXh0U291cmNlXG5cdCkge1xuXHRcdG5vZGUuX3NvdXJjZS5fdW5zdWJzY3JpYmUobm9kZSk7XG5cdH1cblx0ZWZmZWN0Ll9mbiA9IHVuZGVmaW5lZDtcblx0ZWZmZWN0Ll9zb3VyY2VzID0gdW5kZWZpbmVkO1xuXG5cdGNsZWFudXBFZmZlY3QoZWZmZWN0KTtcbn1cblxuZnVuY3Rpb24gZW5kRWZmZWN0KHRoaXM6IEVmZmVjdCwgcHJldkNvbnRleHQ/OiBDb21wdXRlZCB8IEVmZmVjdCkge1xuXHRpZiAoZXZhbENvbnRleHQgIT09IHRoaXMpIHtcblx0XHR0aHJvdyBuZXcgRXJyb3IoXCJPdXQtb2Ytb3JkZXIgZWZmZWN0XCIpO1xuXHR9XG5cdGNsZWFudXBTb3VyY2VzKHRoaXMpO1xuXHRldmFsQ29udGV4dCA9IHByZXZDb250ZXh0O1xuXG5cdHRoaXMuX2ZsYWdzICY9IH5SVU5OSU5HO1xuXHRpZiAodGhpcy5fZmxhZ3MgJiBESVNQT1NFRCkge1xuXHRcdGRpc3Bvc2VFZmZlY3QodGhpcyk7XG5cdH1cblx0ZW5kQmF0Y2goKTtcbn1cblxudHlwZSBFZmZlY3RGbiA9XG5cdHwgKCh0aGlzOiB7IGRpc3Bvc2U6ICgpID0+IHZvaWQgfSkgPT4gdm9pZCB8ICgoKSA9PiB2b2lkKSlcblx0fCAoKCkgPT4gdm9pZCB8ICgoKSA9PiB2b2lkKSk7XG5cbi8vIEF2b2lkIGhhcmQtcmVxdWlyaW5nIHRoZSBFU05leHQuRGlzcG9zYWJsZSBsaWIgaW4gY29uc3VtaW5nIHRzY29uZmlncy5cbi8vIFdoZW4gYFN5bWJvbC5kaXNwb3NlYCBpcyBhdmFpbGFibGUsIHRoaXMgYmVjb21lcyBhIHN5bWJvbC1rZXllZCBkaXNwb3NlciB0eXBlLlxudHlwZSBEaXNwb3NlU3ltYm9sID0gdHlwZW9mIFN5bWJvbCBleHRlbmRzIHsgcmVhZG9ubHkgZGlzcG9zZTogaW5mZXIgVERpc3Bvc2UgfVxuXHQ/IFREaXNwb3NlXG5cdDogbmV2ZXI7XG50eXBlIERpc3Bvc2FibGVMaWtlID0ge1xuXHRbSyBpbiBEaXNwb3NlU3ltYm9sICYgUHJvcGVydHlLZXldOiAoKSA9PiB2b2lkO1xufTtcbnR5cGUgRGlzcG9zZUZuID0gKCgpID0+IHZvaWQpICYgRGlzcG9zYWJsZUxpa2U7XG5cbi8qKlxuICogVGhlIGJhc2UgY2xhc3MgZm9yIHJlYWN0aXZlIGVmZmVjdHMuXG4gKi9cbmRlY2xhcmUgY2xhc3MgRWZmZWN0IHtcblx0X2ZuPzogRWZmZWN0Rm47XG5cdF9jbGVhbnVwPzogKCkgPT4gdm9pZDtcblx0X3NvdXJjZXM/OiBOb2RlO1xuXHRfbmV4dEJhdGNoZWRFZmZlY3Q/OiBFZmZlY3Q7XG5cdF9mbGFnczogbnVtYmVyO1xuXHRfZGVidWdDYWxsYmFjaz86ICgpID0+IHZvaWQ7XG5cdG5hbWU/OiBzdHJpbmc7XG5cblx0Y29uc3RydWN0b3IoZm46IEVmZmVjdEZuLCBvcHRpb25zPzogRWZmZWN0T3B0aW9ucyk7XG5cblx0X2NhbGxiYWNrKCk6IHZvaWQ7XG5cdF9zdGFydCgpOiAoKSA9PiB2b2lkO1xuXHRfbm90aWZ5KCk6IHZvaWQ7XG5cdF9kaXNwb3NlKCk6IHZvaWQ7XG5cdGRpc3Bvc2UoKTogdm9pZDtcbn1cblxuZXhwb3J0IGludGVyZmFjZSBFZmZlY3RPcHRpb25zIHtcblx0bmFtZT86IHN0cmluZztcbn1cblxuLyoqIEBpbnRlcm5hbCAqL1xuZnVuY3Rpb24gRWZmZWN0KHRoaXM6IEVmZmVjdCwgZm46IEVmZmVjdEZuLCBvcHRpb25zPzogRWZmZWN0T3B0aW9ucykge1xuXHR0aGlzLl9mbiA9IGZuO1xuXHR0aGlzLl9jbGVhbnVwID0gdW5kZWZpbmVkO1xuXHR0aGlzLl9zb3VyY2VzID0gdW5kZWZpbmVkO1xuXHR0aGlzLl9uZXh0QmF0Y2hlZEVmZmVjdCA9IHVuZGVmaW5lZDtcblx0dGhpcy5fZmxhZ3MgPSBUUkFDS0lORztcblx0dGhpcy5uYW1lID0gb3B0aW9ucz8ubmFtZTtcblxuXHRpZiAoY2FwdHVyZWRFZmZlY3RzKSB7XG5cdFx0Y2FwdHVyZWRFZmZlY3RzLnB1c2godGhpcyk7XG5cdH1cbn1cblxuRWZmZWN0LnByb3RvdHlwZS5fY2FsbGJhY2sgPSBmdW5jdGlvbiAoKSB7XG5cdGNvbnN0IGZpbmlzaCA9IHRoaXMuX3N0YXJ0KCk7XG5cdHRyeSB7XG5cdFx0aWYgKHRoaXMuX2ZsYWdzICYgRElTUE9TRUQpIHJldHVybjtcblx0XHRpZiAodGhpcy5fZm4gPT09IHVuZGVmaW5lZCkgcmV0dXJuO1xuXG5cdFx0Y29uc3QgY2xlYW51cCA9IHRoaXMuX2ZuKCk7XG5cdFx0aWYgKHR5cGVvZiBjbGVhbnVwID09PSBcImZ1bmN0aW9uXCIpIHtcblx0XHRcdHRoaXMuX2NsZWFudXAgPSBjbGVhbnVwO1xuXHRcdH1cblx0fSBmaW5hbGx5IHtcblx0XHRmaW5pc2goKTtcblx0fVxufTtcblxuRWZmZWN0LnByb3RvdHlwZS5fc3RhcnQgPSBmdW5jdGlvbiAoKSB7XG5cdGlmICh0aGlzLl9mbGFncyAmIFJVTk5JTkcpIHtcblx0XHR0aHJvdyBuZXcgRXJyb3IoXCJDeWNsZSBkZXRlY3RlZFwiKTtcblx0fVxuXHR0aGlzLl9mbGFncyB8PSBSVU5OSU5HO1xuXHR0aGlzLl9mbGFncyAmPSB+RElTUE9TRUQ7XG5cdGNsZWFudXBFZmZlY3QodGhpcyk7XG5cdHByZXBhcmVTb3VyY2VzKHRoaXMpO1xuXG5cdC8qQF9fSU5MSU5FX18qKi8gc3RhcnRCYXRjaCgpO1xuXHRjb25zdCBwcmV2Q29udGV4dCA9IGV2YWxDb250ZXh0O1xuXHRldmFsQ29udGV4dCA9IHRoaXM7XG5cdHJldHVybiBlbmRFZmZlY3QuYmluZCh0aGlzLCBwcmV2Q29udGV4dCk7XG59O1xuXG5FZmZlY3QucHJvdG90eXBlLl9ub3RpZnkgPSBmdW5jdGlvbiAoKSB7XG5cdGlmICghKHRoaXMuX2ZsYWdzICYgTk9USUZJRUQpKSB7XG5cdFx0dGhpcy5fZmxhZ3MgfD0gTk9USUZJRUQ7XG5cdFx0dGhpcy5fbmV4dEJhdGNoZWRFZmZlY3QgPSBiYXRjaGVkRWZmZWN0O1xuXHRcdGJhdGNoZWRFZmZlY3QgPSB0aGlzO1xuXHR9XG59O1xuXG5FZmZlY3QucHJvdG90eXBlLl9kaXNwb3NlID0gZnVuY3Rpb24gKCkge1xuXHR0aGlzLl9mbGFncyB8PSBESVNQT1NFRDtcblxuXHRpZiAoISh0aGlzLl9mbGFncyAmIFJVTk5JTkcpKSB7XG5cdFx0ZGlzcG9zZUVmZmVjdCh0aGlzKTtcblx0fVxufTtcblxuRWZmZWN0LnByb3RvdHlwZS5kaXNwb3NlID0gZnVuY3Rpb24gKCkge1xuXHR0aGlzLl9kaXNwb3NlKCk7XG59O1xuLyoqXG4gKiBDcmVhdGUgYW4gZWZmZWN0IHRvIHJ1biBhcmJpdHJhcnkgY29kZSBpbiByZXNwb25zZSB0byBzaWduYWwgY2hhbmdlcy5cbiAqXG4gKiBBbiBlZmZlY3QgdHJhY2tzIHdoaWNoIHNpZ25hbHMgYXJlIGFjY2Vzc2VkIHdpdGhpbiB0aGUgZ2l2ZW4gY2FsbGJhY2tcbiAqIGZ1bmN0aW9uIGBmbmAsIGFuZCByZS1ydW5zIHRoZSBjYWxsYmFjayB3aGVuIHRob3NlIHNpZ25hbHMgY2hhbmdlLlxuICpcbiAqIFRoZSBjYWxsYmFjayBtYXkgcmV0dXJuIGEgY2xlYW51cCBmdW5jdGlvbi4gVGhlIGNsZWFudXAgZnVuY3Rpb24gZ2V0c1xuICogcnVuIG9uY2UsIGVpdGhlciB3aGVuIHRoZSBjYWxsYmFjayBpcyBuZXh0IGNhbGxlZCBvciB3aGVuIHRoZSBlZmZlY3RcbiAqIGdldHMgZGlzcG9zZWQsIHdoaWNoZXZlciBoYXBwZW5zIGZpcnN0LlxuICpcbiAqIEBwYXJhbSBmbiBUaGUgZWZmZWN0IGNhbGxiYWNrLlxuICogQHJldHVybnMgQSBmdW5jdGlvbiBmb3IgZGlzcG9zaW5nIHRoZSBlZmZlY3QuXG4gKi9cbmZ1bmN0aW9uIGVmZmVjdChmbjogRWZmZWN0Rm4sIG9wdGlvbnM/OiBFZmZlY3RPcHRpb25zKTogRGlzcG9zZUZuIHtcblx0Y29uc3QgZWZmZWN0ID0gbmV3IEVmZmVjdChmbiwgb3B0aW9ucyk7XG5cdHRyeSB7XG5cdFx0ZWZmZWN0Ll9jYWxsYmFjaygpO1xuXHR9IGNhdGNoIChlcnIpIHtcblx0XHRlZmZlY3QuX2Rpc3Bvc2UoKTtcblx0XHR0aHJvdyBlcnI7XG5cdH1cblx0Ly8gUmV0dXJuIGEgYm91bmQgZnVuY3Rpb24gaW5zdGVhZCBvZiBhIHdyYXBwZXIgbGlrZSBgKCkgPT4gZWZmZWN0Ll9kaXNwb3NlKClgLFxuXHQvLyBiZWNhdXNlIGJvdW5kIGZ1bmN0aW9ucyBzZWVtIHRvIGJlIGp1c3QgYXMgZmFzdCBhbmQgdGFrZSB1cCBhIGxvdCBsZXNzIG1lbW9yeS5cblx0Y29uc3QgZGlzcG9zZSA9IGVmZmVjdC5fZGlzcG9zZS5iaW5kKGVmZmVjdCk7XG5cdChkaXNwb3NlIGFzIGFueSlbU3ltYm9sLmRpc3Bvc2VdID0gZGlzcG9zZTtcblx0cmV0dXJuIGRpc3Bvc2UgYXMgRGlzcG9zZUZuO1xufVxuXG4vLyNlbmRyZWdpb24gRWZmZWN0XG5cbi8vI3JlZ2lvbiBBY3Rpb25cblxuZnVuY3Rpb24gYWN0aW9uPFRBcmdzIGV4dGVuZHMgdW5rbm93bltdLCBUUmV0dXJuPihcblx0Zm46ICguLi5hcmdzOiBUQXJncykgPT4gVFJldHVyblxuKTogKC4uLmFyZ3M6IFRBcmdzKSA9PiBUUmV0dXJuIHtcblx0cmV0dXJuIGZ1bmN0aW9uIGFjdGlvbldyYXBwZXIodGhpczogdW5rbm93biwgLi4uYXJnczogVEFyZ3MpIHtcblx0XHRyZXR1cm4gYmF0Y2goKCkgPT4gdW50cmFja2VkKCgpID0+IGZuLmFwcGx5KHRoaXMsIGFyZ3MpKSk7XG5cdH07XG59XG5cbi8vI2VuZHJlZ2lvbiBBY3Rpb25cblxuLy8jcmVnaW9uIGNyZWF0ZU1vZGVsXG5cbi8qKiBNb2RlbHMgc2hvdWxkIG9ubHkgY29udGFpbiBzaWduYWxzLCBhY3Rpb25zLCBhbmQgbmVzdGVkIG9iamVjdHMgY29udGFpbmluZyBvbmx5IHNpZ25hbHMgYW5kIGFjdGlvbnMuICovXG50eXBlIFZhbGlkYXRlTW9kZWw8VE1vZGVsPiA9IHtcblx0W0tleSBpbiBrZXlvZiBUTW9kZWxdOiBUTW9kZWxbS2V5XSBleHRlbmRzIFJlYWRvbmx5U2lnbmFsPHVua25vd24+XG5cdFx0PyBUTW9kZWxbS2V5XVxuXHRcdDogVE1vZGVsW0tleV0gZXh0ZW5kcyAoLi4uYXJnczogYW55W10pID0+IGFueVxuXHRcdFx0PyBUTW9kZWxbS2V5XVxuXHRcdFx0OiBUTW9kZWxbS2V5XSBleHRlbmRzIG9iamVjdFxuXHRcdFx0XHQ/IFZhbGlkYXRlTW9kZWw8VE1vZGVsW0tleV0+XG5cdFx0XHRcdDogYFByb3BlcnR5ICR7S2V5IGV4dGVuZHMgc3RyaW5nID8gYCcke0tleX0nIGAgOiBcIlwifWlzIG5vdCBhIFNpZ25hbCwgQWN0aW9uLCBvciBhbiBvYmplY3QgdGhhdCBjb250YWlucyBvbmx5IFNpZ25hbHMgYW5kIEFjdGlvbnMuYDtcbn07XG5cbmV4cG9ydCB0eXBlIE1vZGVsPFRNb2RlbD4gPSBWYWxpZGF0ZU1vZGVsPFRNb2RlbD4gJiBEaXNwb3NhYmxlTGlrZTtcblxuZXhwb3J0IHR5cGUgTW9kZWxGYWN0b3J5PFRNb2RlbCwgVEZhY3RvcnlBcmdzIGV4dGVuZHMgYW55W10gPSBbXT4gPSAoXG5cdC4uLmFyZ3M6IFRGYWN0b3J5QXJnc1xuKSA9PiBWYWxpZGF0ZU1vZGVsPFRNb2RlbD47XG5leHBvcnQgdHlwZSBNb2RlbENvbnN0cnVjdG9yPFRNb2RlbCwgVEZhY3RvcnlBcmdzIGV4dGVuZHMgYW55W10gPSBbXT4gPSBuZXcgKFxuXHQuLi5hcmdzOiBURmFjdG9yeUFyZ3NcbikgPT4gTW9kZWw8VE1vZGVsPjtcblxuLyoqXG4gKiBUaGUgcHVibGljIHR5cGVzIGZvciBNb2RlbENvbnN0cnVjdG9yIHJlcXVpcmUgdXNpbmcgYG5ld2AgdG8gaGVscFxuICogZGlzYW1iaWd1YXRlIHRoZSBmdW5jdGlvbiBwYXNzZWQgaW50byBgY3JlYXRlTW9kZWxgIGFuZCB0aGUgcmV0dXJuZWRcbiAqIGNvbnN0cnVjdG9yIGZ1bmN0aW9uLiBJdCBpcyBlYXNpZXIgdG8gc2F5IHRoYXQgYGNyZWF0ZU1vZGVsYCBhY2NlcHRzXG4gKiBhIGZhY3RvcnkgYW5kIHJldHVybnMgYSBjbGFzcywgdGhlbiB0byBzYXkgaXQgYWNjZXB0cyBhIGZhY3RvcnkgYW5kXG4gKiByZXR1cm5zIGEgZmFjdG9yeS4gSW4gb3RoZXIgd29yZHMsIHRoaXMgZXhhbXBsZTpcbiAqXG4gKiBgYGB0c1xuICogY29uc3QgUGVyc29uTW9kZWwgPSBjcmVhdGVNb2RlbCgobmFtZTogc3RyaW5nKSA9PiAoeyAuLi4gfSkpO1xuICogY29uc3QgcGVyc29uID0gbmV3IFBlcnNvbk1vZGVsKFwiSm9oblwiKTtcbiAqIGBgYFxuICpcbiAqIGlzIGVhc2llciB0byB1bmRlcnN0YW5kIHRoYW4gdGhpcyBleGFtcGxlOlxuICpcbiAqIGBgYHRzXG4gKiBjb25zdCBjcmVhdGVQZXJzb24gPSBjcmVhdGVNb2RlbCgobmFtZTogc3RyaW5nKSA9PiAoeyAuLi4gfSkpO1xuICogY29uc3QgcGVyc29uID0gY3JlYXRlUGVyc29uKFwiSm9oblwiKTtcbiAqIGBgYFxuICpcbiAqIEhvd2V2ZXIsIGludGVybmFsbHkgd2UgaW1wbGVtZW50IGBjcmVhdGVNb2RlbGAgdG8gcmV0dXJuIGEgZnVuY3Rpb25cbiAqIHRoYXQgY2FuIGJlIGNhbGxlZCB3aXRob3V0IGBuZXdgIGZvciBzaW1wbGljaXR5LiBUbyBicmlkZ2UgdGhlIGdhcFxuICogYmV0d2VlbiB0aGUgcHVibGljIHR5cGVzIGFuZCB0aGUgaW50ZXJuYWwgaW1wbGVtZW50YXRpb24sIHdlIGRlZmluZVxuICogdGhpcyBpbnRlcm5hbCBpbnRlcmZhY2UgdGhhdCBleHRlbmRzIHRoZSBwdWJsaWMgaW50ZXJmYWNlIGJ1dCBhbHNvXG4gKiBhbGxvd3MgY2FsbGluZyB3aXRob3V0IGBuZXdgLlxuICpcbiAqIFRoaXMgcGF0dGVybiBpcyB1c2VkIGJ5IHRoZSBQcmVhY3QgJiBSZWFjdCBhZGFwdGVycyB0byBtYWtlIGluc3RhbnRpYXRpbmdcbiAqIGEgbW9kZWwgb3IgYSBmdW5jdGlvbiB0aGF0IHJldHVybnMgYSBtb2RlbCBlYXNpZXIuXG4gKlxuICogQGludGVybmFsXG4gKi9cbmludGVyZmFjZSBJbnRlcm5hbE1vZGVsQ29uc3RydWN0b3I8XG5cdFRNb2RlbCxcblx0VEZhY3RvcnlBcmdzIGV4dGVuZHMgYW55W10sXG4+IGV4dGVuZHMgTW9kZWxDb25zdHJ1Y3RvcjxUTW9kZWwsIFRGYWN0b3J5QXJncz4ge1xuXHQoLi4uYXJnczogVEZhY3RvcnlBcmdzKTogTW9kZWw8VE1vZGVsPjtcbn1cblxuZnVuY3Rpb24gc3RhcnRDYXB0dXJpbmdFZmZlY3RzKCk6ICgpID0+IEVmZmVjdFtdIHwgdW5kZWZpbmVkIHtcblx0bGV0IHByZXZDYXB0dXJlZEVmZmVjdHMgPSBjYXB0dXJlZEVmZmVjdHM7XG5cdC8vIEFsd2F5cyBlc3RhYmxpc2ggYSBmcmVzaCBjYXB0dXJlIHNjb3BlLCBldmVuIHdoZW4gYHVudHJhY2tlZCgpYCBoYXNcblx0Ly8gdGVtcG9yYXJpbHkgY2xlYXJlZCB0aGUgcGFyZW50IHNjb3BlLiBUaGlzIGxldHMgbmVzdGVkIG1vZGVscyBvd24gdGhlaXJcblx0Ly8gZWZmZWN0cyB3aXRob3V0IHByb21vdGluZyB0aGVtIHRvIGEgc3VwcHJlc3NlZCBvdXRlciBzY29wZS5cblx0Y2FwdHVyZWRFZmZlY3RzID0gW107XG5cblx0cmV0dXJuIGZ1bmN0aW9uIHN0b3BDYXB0dXJpbmdFZmZlY3RzKCkge1xuXHRcdGxldCBtb2RlbEVmZmVjdHMgPSBjYXB0dXJlZEVmZmVjdHM7XG5cdFx0aWYgKGNhcHR1cmVkRWZmZWN0cyAmJiBwcmV2Q2FwdHVyZWRFZmZlY3RzKSB7XG5cdFx0XHRwcmV2Q2FwdHVyZWRFZmZlY3RzID0gcHJldkNhcHR1cmVkRWZmZWN0cy5jb25jYXQoY2FwdHVyZWRFZmZlY3RzKTtcblx0XHR9XG5cblx0XHRjYXB0dXJlZEVmZmVjdHMgPSBwcmV2Q2FwdHVyZWRFZmZlY3RzO1xuXG5cdFx0cmV0dXJuIG1vZGVsRWZmZWN0cztcblx0fTtcbn1cblxuY29uc3Qgd3JhcEluQWN0aW9uID0gKHZhbHVlOiBSZWNvcmQ8c3RyaW5nLCB1bmtub3duPikgPT4ge1xuXHRmb3IgKGNvbnN0IGtleSBpbiB2YWx1ZSkge1xuXHRcdGNvbnN0IHZhbCA9IHZhbHVlW2tleV07XG5cdFx0aWYgKHR5cGVvZiB2YWwgPT09IFwiZnVuY3Rpb25cIikge1xuXHRcdFx0dmFsdWVba2V5XSA9IGFjdGlvbih2YWwgYXMgKC4uLmFyZ3M6IHVua25vd25bXSkgPT4gdW5rbm93bik7XG5cdFx0fSBlbHNlIGlmICh0eXBlb2YgdmFsID09PSBcIm9iamVjdFwiICYmIHZhbCAhPT0gbnVsbCAmJiAhKFwiYnJhbmRcIiBpbiB2YWwpKSB7XG5cdFx0XHQvLyBSZWN1cnNpdmVseSB3cmFwIG5lc3RlZCBvYmplY3QgcHJvcGVydGllcyBpbiBhY3Rpb25zLiBUaGlzIGFsbG93cyB1c2VycyB0byB3cml0ZVxuXHRcdFx0Ly8gbmVzdGVkIG1vZGVscyB3aXRob3V0IHdvcnJ5aW5nIGFib3V0IHdyYXBwaW5nIHRoZWlyIGZ1bmN0aW9ucyBpbiBgYWN0aW9uYC5cblx0XHRcdHdyYXBJbkFjdGlvbih2YWwgYXMgUmVjb3JkPHN0cmluZywgdW5rbm93bj4pO1xuXHRcdH1cblx0fVxufTtcblxuZnVuY3Rpb24gY3JlYXRlTW9kZWw8VE1vZGVsLCBURmFjdG9yeUFyZ3MgZXh0ZW5kcyBhbnlbXSA9IFtdPihcblx0bW9kZWxGYWN0b3J5OiBNb2RlbEZhY3Rvcnk8VE1vZGVsLCBURmFjdG9yeUFyZ3M+XG4pOiBNb2RlbENvbnN0cnVjdG9yPFRNb2RlbCwgVEZhY3RvcnlBcmdzPiB7XG5cdHJldHVybiBmdW5jdGlvbiBTaWduYWxNb2RlbCguLi5hcmdzOiBURmFjdG9yeUFyZ3MpOiBNb2RlbDxUTW9kZWw+IHtcblx0XHRsZXQgbW9kZWxFZmZlY3RzOiBFZmZlY3RbXSB8IHVuZGVmaW5lZDtcblx0XHRsZXQgbW9kZWw6IE1vZGVsPFRNb2RlbD47XG5cblx0XHRjb25zdCBzdG9wQ2FwdHVyaW5nRWZmZWN0cyA9IHN0YXJ0Q2FwdHVyaW5nRWZmZWN0cygpO1xuXHRcdHRyeSB7XG5cdFx0XHRtb2RlbCA9IG1vZGVsRmFjdG9yeSguLi5hcmdzKSBhcyBNb2RlbDxUTW9kZWw+O1xuXHRcdH0gY2F0Y2ggKGVycikge1xuXHRcdFx0Ly8gRHJvcCBhbnkgY2FwdHVyZWQgZWZmZWN0cyBvbiBlcnJvci4gRXJyb3JzIGZyb20gbmVzdGVkIG1vZGVscyB3aWxsIGJ1YmJsZVxuXHRcdFx0Ly8gdXAgaGVyZSBhbmQgcmVjdXJzaXZlbHkgcmVzZXQgYGNhcHR1cmVkRWZmZWN0c2AgdG8gYHVuZGVmaW5lZGAgcHJldmVudGluZ1xuXHRcdFx0Ly8gYW55IGNhcHR1cmVkIGVmZmVjdHMgZnJvbSBsZWFraW5nXG5cdFx0XHRjYXB0dXJlZEVmZmVjdHMgPSB1bmRlZmluZWQ7XG5cdFx0XHR0aHJvdyBlcnI7XG5cdFx0fSBmaW5hbGx5IHtcblx0XHRcdG1vZGVsRWZmZWN0cyA9IHN0b3BDYXB0dXJpbmdFZmZlY3RzKCk7XG5cdFx0fVxuXG5cdFx0d3JhcEluQWN0aW9uKG1vZGVsKTtcblxuXHRcdG1vZGVsW1N5bWJvbC5kaXNwb3NlXSA9IGFjdGlvbihmdW5jdGlvbiBkaXNwb3NlTW9kZWwoKSB7XG5cdFx0XHRpZiAobW9kZWxFZmZlY3RzKSB7XG5cdFx0XHRcdGZvciAobGV0IGkgPSAwOyBpIDwgbW9kZWxFZmZlY3RzLmxlbmd0aDsgaSsrKSB7XG5cdFx0XHRcdFx0bW9kZWxFZmZlY3RzW2ldLmRpc3Bvc2UoKTtcblx0XHRcdFx0fVxuXHRcdFx0fVxuXG5cdFx0XHRtb2RlbEVmZmVjdHMgPSB1bmRlZmluZWQ7XG5cdFx0fSk7XG5cblx0XHRyZXR1cm4gbW9kZWw7XG5cdH0gYXMgSW50ZXJuYWxNb2RlbENvbnN0cnVjdG9yPFRNb2RlbCwgVEZhY3RvcnlBcmdzPjtcbn1cblxuLy8jZW5kcmVnaW9uIGNyZWF0ZU1vZGVsXG5cbmV4cG9ydCB7XG5cdGNvbXB1dGVkLFxuXHRlZmZlY3QsXG5cdGJhdGNoLFxuXHR1bnRyYWNrZWQsXG5cdGFjdGlvbixcblx0Y3JlYXRlTW9kZWwsXG5cdFNpZ25hbCxcblx0UmVhZG9ubHlTaWduYWwsXG5cdEVmZmVjdCxcblx0Q29tcHV0ZWQsXG59O1xuIl19