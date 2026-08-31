/* The shared 'Appearance' prop group (APPEARANCE.md): 13 root-safe CSS props every visual tag
   answers, injected at registration and applied once after create. Absent means platform styling —
   there are no defaults, and nothing is ever serialized for an unassigned prop. */
import {Signal} from '../core/signals.js';
import type {Scope} from '../core/scope.js';
import {Control} from '../core/component.js';
import type {SpecPropMeta} from './registry.js';

export const APPEARANCE_CATEGORY = 'Appearance';

export const APPEARANCE_PROPS: SpecPropMeta[] = ([
  {name: 'color', type: 'string', inputType: 'Color', description: 'Text color; inherits down the subtree.'},
  {name: 'backgroundColor', type: 'string', inputType: 'Color', description: 'Background color.'},
  {name: 'font', type: 'string', inputType: 'Font', description: 'CSS font shorthand ("bold 14px Roboto").'},
  {name: 'opacity', type: 'double', min: 0, max: 1, step: 0.05, description: '0 transparent to 1 opaque.'},
  {name: 'cursor', type: 'string', nullable: true,
    choices: ['default', 'pointer', 'not-allowed', 'wait', 'help', 'text', 'move', 'grab'],
    description: 'Mouse cursor over the node.'},
  {name: 'textAlign', type: 'string', nullable: true, choices: ['left', 'center', 'right'],
    description: 'Horizontal text alignment.'},
  {name: 'padding', type: 'string', description: 'CSS padding ("8px 12px").'},
  {name: 'margin', type: 'string', description: 'CSS margin.'},
  {name: 'border', type: 'string', description: 'CSS border ("1px solid var(--dg-border)").'},
  {name: 'borderRadius', type: 'string', description: 'Corner radius ("4px").'},
  {name: 'width', type: 'string', description: 'CSS width ("200px", "50%").'},
  {name: 'height', type: 'string', description: 'CSS height.'},
  {name: 'pointerEvents', type: 'string', nullable: true, choices: ['auto', 'none'],
    description: '"none" makes the node click-through.'},
] as SpecPropMeta[]).map((p) => ({...p, category: APPEARANCE_CATEGORY, bindable: true}));

/** The color-token palette the designer's ColorInput offers, semantic entries first — full
 * `var(--dg-…)` strings, so the token lint validates every entry against `css/tokens.css`. */
export const DESIGN_TOKENS: string[] = [
  'var(--dg-text-color)', 'var(--dg-text-color-light)', 'var(--dg-bg)', 'var(--dg-bg-secondary)',
  'var(--dg-accent)', 'var(--dg-success)', 'var(--dg-failure)', 'var(--dg-error)',
  'var(--dg-white)', 'var(--dg-grey-1)', 'var(--dg-grey-2)', 'var(--dg-grey-3)',
  'var(--dg-grey-4)', 'var(--dg-grey-5)', 'var(--dg-grey-6)',
  'var(--dg-steel-1)', 'var(--dg-steel-2)', 'var(--dg-steel-3)', 'var(--dg-steel-4)', 'var(--dg-steel-5)',
  'var(--dg-blue-1)', 'var(--dg-blue-2)',
  'var(--dg-green-1)', 'var(--dg-green-2)', 'var(--dg-green-3)',
  'var(--dg-orange-1)', 'var(--dg-orange-2)', 'var(--dg-orange-3)',
  'var(--dg-red-1)', 'var(--dg-red-3)', 'var(--dg-red-4)',
  'var(--dg-warm-grey-1)', 'var(--dg-warm-grey-2)', 'var(--dg-warm-grey-3)', 'var(--dg-warm-grey-4)',
  'var(--dg-beige-1)', 'var(--dg-beige-2)',
];

/** The shared list minus name collisions — a component's own same-named prop wins. */
export function appearanceFor(props: SpecPropMeta[]): SpecPropMeta[] {
  return APPEARANCE_PROPS.filter((a) => !props.some((p) => p.name === a.name));
}

/** Writes each applicable prop present in `props` onto the target's root style: a literal once, a
 * bound signal through an effect on the control's own scope — it dies with the node. For a plain
 * element a bound signal needs `scope` (the node's mount scope); without one it is skipped, as the
 * registry's element path always was. */
export function applyAppearance(target: Control | HTMLElement, props: Record<string, unknown>,
  applicable: SpecPropMeta[], scope?: Scope): void {
  const control = Control.is(target) ? target : undefined;
  const style = (control?.root ?? target as HTMLElement).style;
  for (const prop of applicable) {
    const name = prop.name;
    const value = props[name];
    if (value === undefined)
      continue;
    if (value instanceof Signal) {
      if (control === undefined) {
        if (scope !== undefined)
          scope.effect(() => write((target as HTMLElement).style, name, (value as Signal<unknown>).value));
        continue;
      }
      /* The bound signal becomes a same-named member so `bindStep(name)` answers it: a deferred
         bind's `link` finds the proxy there and wires the source through it, and `bindProps()`
         advertises exactly the bound appearance props. The `in` guard never clobbers an existing
         member, even a null-valued one. */
      if (!(name in control))
        (control as unknown as Record<string, unknown>)[name] = value;
      // root read at write time: Control has a root setter, and a captured style would go stale
      control.effect(() => write(control.root.style, name, (value as Signal<unknown>).value));
    } else
      write(style, name, value);
  }
}

function write(style: CSSStyleDeclaration, name: string, value: unknown): void {
  (style as unknown as Record<string, unknown>)[name] = value == null || value === '' ? '' : String(value);
}
