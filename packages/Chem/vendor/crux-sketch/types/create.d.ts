import type { CreateSketcherOptions, SketcherApi } from './core/types.js';
/** What `createSketcher` resolves to: the sketcher's API. */
export type SketcherHandle = SketcherApi;
/**
 * Creates a sketcher inside `element` and resolves with its handle once it is ready: the engine is
 * loaded and the initial value, if any, is drawn. Showing the initial value emits no change event;
 * a value the sketcher cannot read leaves it empty, with `lastError` set. If the engine cannot be
 * loaded, it rejects, and leaves nothing in `element`.
 */
export declare function createSketcher(element: Element, options?: CreateSketcherOptions): Promise<SketcherHandle>;
