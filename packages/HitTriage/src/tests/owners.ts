import {Scope} from '@datagrok-libraries/u2';
import type {Component} from '@datagrok-libraries/u2';

/** What a test builds outside any view: [build]'s component on a root of its own (`Scope.root(label)`), ended with
 * it. */
export function onRoot<C extends Component>(label: string, build: () => C): C {
  const owner = Scope.root(label);
  const built = owner.build(build);
  built.own(() => owner.dispose());
  return built;
}
