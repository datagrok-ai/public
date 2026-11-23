/** Preload (NODE_OPTIONS=--import) that makes a loaded machine on purpose: a timer of 15-25 ms (the 20 ms read
   timeout the binding-dialog tests set) fires 150 ms late. A test that sleeps a fixed time past it always loses;
   one that waits for the condition still passes. Loaded only by the EQ-032 regression command, never by the suite. */
const realSetTimeout = globalThis.setTimeout;
globalThis.setTimeout = function(fn, ms, ...args) {
  return realSetTimeout.call(this, fn, ms >= 15 && ms <= 25 ? ms + 150 : ms, ...args);
};
