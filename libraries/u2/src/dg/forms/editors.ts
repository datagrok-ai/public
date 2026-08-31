/* Editor selection for schema-driven forms (OBJECTS.md): rules registered here are consulted
   by ObjectForm before `editors: 'auto'` and the generated editors, so semType-specific inputs
   (molecule, sequence…) plug in without the generator learning any domain. Structural on
   purpose — platform-bound editors register themselves from their own modules. */
import {Input, InputOptions} from '../../core/input-base.js';
import type {IProperty} from '../../core/property-like.js';

export interface EditorRule {
  match(prop: IProperty): boolean;
  /** `options` carries the generator's per-field bag (label, name, tooltip, overrides). The
   * form owns the created input and wires it by the property's declared type, so the editor's
   * value must be that type (a Molecule editor is an `Input<string>`). */
  create(prop: IProperty, options: InputOptions<any>): Input<any>;
}

export class Editors {
  private static _rules: EditorRule[] = [];

  /** First matching rule wins, in registration order. Returns the unregister function. */
  static register(rule: EditorRule): () => void {
    Editors._rules.push(rule);
    return () => {
      const i = Editors._rules.indexOf(rule);
      if (i >= 0)
        Editors._rules.splice(i, 1);
    };
  }

  static resolve(prop: IProperty, options: InputOptions<any>): Input<any> | null {
    for (const rule of Editors._rules) {
      if (rule.match(prop))
        return rule.create(prop, options);
    }
    return null;
  }
}
