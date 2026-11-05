// any control loaded installs the page's one u2 stylesheet
import '../styles.js';
import {Component as EngineComponent, Control as EngineControl, component as engineComponent, version}
  from 'datagrok-api/u2core';
import type {OwnerVerbs, Scope} from 'datagrok-api/u2core';

export type {PropertyChange, ComponentState, OwnerVerbs, BuiltRoot, BuiltHandle,
  ControlBuildOptions} from 'datagrok-api/u2core';

declare module 'datagrok-api/u2core' {
  // the engine puts the owner verbs on every component and u2 declares them: `dialog.build(fn)` type-checks here,
  // while a platform class with a `build` of its own (a function view), compiled without u2, keeps compiling
  interface Component extends OwnerVerbs {}
}

// In a plugin bundle the engine is the platform's DG.U2, which an older platform serves without what this u2 calls.
// The package still loads there; what builds a u2 component refuses, naming the fix.
const engine: string | undefined = version;
const [major, minor] = (engine ?? '').split('.').map(Number);
const refusal = major > 1 || major === 1 && minor >= 2 ? undefined :
  `u2 needs the u2 engine 1.2+; this platform has ${engine ?? 'an older one'}: update the platform`;
const refuse = (): never => {
  throw new Error(refusal);
};

export const Component: typeof EngineComponent = refusal === undefined ? EngineComponent :
  class extends EngineComponent {
    constructor(scope?: Scope) {
      refuse();
      super(scope);
    }
  };
export type Component = EngineComponent;

export const Control: typeof EngineControl = refusal === undefined ? EngineControl :
  class extends EngineControl {
    constructor(root?: HTMLElement, scope?: Scope) {
      refuse();
      super(root, scope);
    }

    static override build(): never {
      return refuse();
    }
  };
export type Control = EngineControl;

export const component: typeof engineComponent = refusal === undefined ? engineComponent : refuse;
