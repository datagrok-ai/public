/* What a domain page is MADE OF, as two functions over a source (Astra C1): the phase-4 designer
   replaces a builder, never a source and never a session — data ownership and the save boundary stay
   with `DomainApp`, which is the whole point of extracting only the assembly. */
import {Scope} from '../../core/scope.js';
import type {Control} from '../../core/component.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {domains} from './index.js';
import {DomainList} from './list.js';
import type {DomainListMode} from './list.js';
import {DomainForm} from './form.js';
import type {DomainChildrenOptions} from './children.js';

export interface BuildListOptions {
  mode?: DomainListMode;
}

export interface BuildEntityOptions {
  /** The form's columns, in this order. */
  include?: string[];
  /** What the form says when the row was not found — the app words it with the address. */
  empty?: string;
  children?: boolean | DomainChildrenOptions;
  history?: boolean;
}

/** The list page's collection over `source`, built in `scope`. */
export function buildList(source: DomainSource, scope: Scope, options: BuildListOptions): DomainList {
  return Scope.runWith(scope, () => new DomainList(source, {mode: options.mode}));
}

/** The entity page over `source`: the form, and the panes under it in their fixed order —
 * the child collections, then the history. `panes` is empty where the options turn both off. */
export function buildEntity(source: DomainSource, scope: Scope, options: BuildEntityOptions):
  {form: DomainForm, panes: Control[]} {
  return Scope.runWith(scope, () => {
    const {children, history} = options;
    return {
      form: new DomainForm(source, {system: 'footer', include: options.include, empty: options.empty}),
      panes: [
        ...(children === false ? [] : [domains.children(source, children === true ? undefined : children)]),
        ...(history === false ? [] : [domains.history(source)]),
      ],
    };
  });
}
