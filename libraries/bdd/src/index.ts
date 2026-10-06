/* What project bindings import from `@datagrok-libraries/bdd`: the registration API and the
   phrase model. Runtime helpers live in `@datagrok-libraries/bdd/runtime`. */
export * from './registry.js';
export * from './nouns.js';
export type {ElementRef} from './runtime/args.js';
export {BASE_TIERS, CONFIG_FILE, PACKAGE_NAME} from './project.js';
export type {Project, ProjectConfig} from './project.js';
