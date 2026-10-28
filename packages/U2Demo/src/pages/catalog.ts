import {
  Control, Registry, SPEC_SCHEMA, SpecContext, Section,
  divV, divH, span, registerAll, renderSpec,
} from '@datagrok-libraries/u2';
import type {ComponentMeta, SpecNode} from '@datagrok-libraries/u2';
import {readout} from './common';
import states from './catalog-states.json';

/** The catalog's states — shared with the generated BDD feature (bdd/tools/gen-catalog.mjs). A state
 * applies to a tag when the tag declares every prop the state `requires`. */
interface CatalogStates {
  order: string[];
  states: Record<string, Record<string, unknown>>;
  requires: Record<string, string[]>;
  /** Tags a state's props cannot put in it (the manifest's `UNREACHED`): a checkbox's `false` is a value. */
  unreached: Record<string, string[]>;
  /** States whose verdict waits for the user (a required field left empty) until the input is touched. */
  touch: string[];
  claims: Record<string, string>;
  layoutOnly: string[];
}

const CATALOG: CatalogStates = states;

/** Every visual registry tag in each state its props reach, rendered from its registry example the way
 * a spec renders it. A private registry holds exactly the tags `registerAll` registers (the manifest's
 * list), whatever else the session has registered globally. */
export function catalogPage(): HTMLElement {
  const page = Control.build(() => {
    const reg = new Registry();
    registerAll(reg);
    let samples = 0;
    let failed = 0;
    const sections: Section[] = [];
    for (const meta of reg.metas()) {
      if (meta.visual === false || CATALOG.layoutOnly.includes(meta.tag))
        continue;
      const short = meta.tag.replace(/^u2-/, '');
      const declared = new Set(meta.props.map((p) => p.name));
      const cells = CATALOG.order
        .filter((state) => (CATALOG.requires[state] ?? []).every((prop) => declared.has(prop)) &&
          !(CATALOG.unreached[state] ?? []).includes(meta.tag))
        .map((state) => {
          const name = `${short}-${state}`;
          let body = sample(reg, meta, name, state);
          const wide = meta.tag === 'u2-splitter' && typeof body !== 'string';
          if (typeof body === 'string') {
            failed++;
            body = span(`failed: ${body}`, 'u2demo-error');
            body.dataset.u2CatalogError = name;
          }
          return divV([span(state, 'u2demo-catalog-state'), body],
            wide ? 'u2demo-catalog-cell u2demo-catalog-cell-wide' : 'u2demo-catalog-cell');
        });
      samples += cells.length;
      sections.push(new Section({title: meta.label ?? short}).add(divH(cells, 'u2demo-catalog-cells')));
    }
    return divV([
      span('Every registered u2 component in each state its props can put it in, rendered from the ' +
        'registry example by the spec renderer. The BDD suite walks this page; a sample that cannot be ' +
        'built says why in red.', 'u2demo-hint'),
      divH([readout('samples', String(samples)), readout('failed', String(failed))], 'u2demo-row'),
      ...sections.map((s) => s.root),
    ], 'u2demo-page u2demo-wide');
  });
  page.name = 'catalog';
  return page.root;
}

/** The sample's root, or the reason it could not be built: the spec renderer reports a node it could
 * not build as an error placeholder rather than a throw. */
function sample(reg: Registry, meta: ComponentMeta, name: string, state: string): HTMLElement | string {
  const example = meta.example as SpecNode;
  // a recipe's null unsets the example's prop (`value: null` — left empty): the renderer type-checks
  // a present value and rejects null
  const props = Object.fromEntries(Object.entries({...example.props, ...CATALOG.states[state]})
    .filter(([, value]) => value !== null));
  try {
    const instance = renderSpec({$schema: SPEC_SCHEMA, root: {...example, name, props}}, new SpecContext(), reg);
    const placeholder = instance.root.querySelector('.u2-spec-error');
    if (placeholder) {
      instance.dispose();
      return placeholder.textContent ?? 'not built';
    }
    if (CATALOG.touch.includes(state))
      (instance.node(name) as {touch?: () => void} | undefined)?.touch?.();
    return instance.root;
  } catch (e) {
    return e instanceof Error ? e.message : String(e);
  }
}
