// Namespace re-exports to avoid name collisions between submodules.
// Usage:
//   import {singleObjective, multiObjective} from './optimization';
//   const nm = new singleObjective.NelderMead();
//   const moead = new multiObjective.Moead(...);
//
// Direct imports from submodules still work:
//   import {NelderMead} from './optimization/single-objective';

export * as singleObjective from './single-objective';
export * as multiObjective from './multi-objectives';
