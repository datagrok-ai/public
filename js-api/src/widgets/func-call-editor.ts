/**
 * FuncCallEditor abstract base class.
 * @module widgets/func-call-editor
 */

import {Observable} from "rxjs";
import {Widget} from "./base";
import type {InputBase} from "./inputs-base";

/**
 * Base class for widgets or views that serve as editors for `FuncCall`. Extend it and use it in editor functions.
 * Editor functions should return an implementation of this class for the platform to handle validation correctly.
 * An editor function can be attached to another function using the `editor` tag: `editor: Plugin:EditorFuncName`.
 */
export abstract class FuncCallEditor extends Widget {
  abstract get isValid(): boolean;
  abstract get onInputChanged(): Observable<any>;

  inputFor?(propertyName: string): InputBase;
}
