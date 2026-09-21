import {FilterBuilder} from '../../components/filter/filter-builder.js';
import {Editors} from '../forms/editors.js';
// side-effect only: loading the schema-driven router is what sets `Editors.byHint`
import '../forms/object-form.js';

export {FilterSchemas} from './schemas.js';
export type {DataFrameFilterSchema} from './schemas.js';
export {frameLike, toBitSet} from './bitset.js';

FilterBuilder.defaultEditors = (prop, options) => Editors.resolve(prop, options);
