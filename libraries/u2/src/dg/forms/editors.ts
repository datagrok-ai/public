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

/** How a non-empty value of a matching property is shown instead of its text: a readonly field, a card, a cell. */
export interface DisplayRule {
  match(prop: IProperty): boolean;
  render(value: unknown, size: DisplaySize, prop: IProperty): HTMLElement;
}

export interface DisplaySize {
  width: number;
  height: number;
}

export class Editors {
  private static _rules: EditorRule[] = [];
  private static _displays: DisplayRule[] = [];
  /** The property's own `inputType`/`editor` hints — the fallback after the rules, set by the
   * schema-driven router so every consumer resolves in the same order. */
  static byHint: ((prop: IProperty, options: InputOptions<any>) => Input<any> | null) | null = null;

  /** First matching rule wins, in registration order. Returns the unregister function. */
  static register(rule: EditorRule): () => void {
    Editors._rules.push(rule);
    return () => {
      const i = Editors._rules.indexOf(rule);
      if (i >= 0)
        Editors._rules.splice(i, 1);
    };
  }

  /** The first matching rule, else what the property's hints ask for, else null. */
  static resolve(prop: IProperty, options: InputOptions<any>): Input<any> | null {
    for (const rule of Editors._rules) {
      if (rule.match(prop))
        return rule.create(prop, options);
    }
    return Editors.byHint?.(prop, options) ?? null;
  }

  /** First matching rule wins, in registration order. Returns the unregister function. */
  static registerDisplay(rule: DisplayRule): () => void {
    Editors._displays.push(rule);
    return () => {
      const i = Editors._displays.indexOf(rule);
      if (i >= 0)
        Editors._displays.splice(i, 1);
    };
  }

  static displays(prop: IProperty): boolean {
    return Editors._displays.some((r) => r.match(prop));
  }

  /** The first matching display rule's rendering of [value], else null. */
  static display(prop: IProperty, value: unknown, size: DisplaySize = {width: 250, height: 160}): HTMLElement | null {
    const rule = Editors._displays.find((r) => r.match(prop));
    return rule === undefined ? null : rule.render(value, size, prop);
  }
}
