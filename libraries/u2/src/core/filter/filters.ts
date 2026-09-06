/* The `Filters` namespace: every module's exports in one flat name space — a name two modules
   share is a compile error here. */
export * from './model.js';
export * from './kinds.js';
export * from './operators.js';
export * from './schema.js';
export * from './validate.js';
export * from './grammar.js';
export * from './domain-tree.js';
export * from './mask.js';
export * from './url.js';
export {resolveSpan, spanOf} from '../span.js';
