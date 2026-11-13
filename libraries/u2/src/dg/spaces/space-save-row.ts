import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {signal, computed} from '../../core/signals.js';
import {component} from '../../core/component.js';
import {div, divH} from '../../core/elements.js';
import {ChoiceChip} from '../../components/inputs/choice-chip.js';
import {statusLine} from '../../components/display/badge.js';
import {SpaceInput} from '../inputs/space-input.js';
import {spaceHasName, rememberSpace} from './spaces.js';

export const SPACE_TOOLTIP = 'Empty: saved to your own projects, visible only to you until you share it.\n' +
  'Pick a space to save it there; it gets that space\'s permissions.';

/** Where OK sends the project: a space id (null for the user's own projects) and whether it is a link. */
export interface SpaceSaveTarget {
  spaceId: string | null;
  link: boolean;
  /** False when the project stays where it is. */
  changed: boolean;
}

/** The Space row of the Save project dialog: the dialog reports the name, holds OK while `problem` is set, asks
 * `resolveTarget()` on OK, then moves or links the project itself. */
export interface SpaceSaveRow extends DG.Widget {
  /** Called whenever `problem` or the target changes. */
  onChanged: () => void;
  /** Why OK cannot save into the picked space now; null when it can. */
  readonly problem: string | null;
  setName(name: string): void;
  /** Where the project goes. */
  resolveTarget(): Promise<SpaceSaveTarget>;
  target: DG.Project | null;
  mode: 'move' | 'link';
}

const same = (a: DG.Project | null, b: DG.Project | null): boolean =>
  a === null || b === null ? a === b : a.id === b.id;

/** `home` is the space a saved project lives in; none for a new project, which always moves. */
export function spaceSaveRow(home: DG.Project | null, name: string): SpaceSaveRow {
  const widget = new DG.Widget(div()) as SpaceSaveRow;
  widget.onChanged = () => {};
  const row = widget.scope.build(() => component((self) => {
    const projectName = signal(name);
    const problem = signal<string | null>(null);
    const input = new SpaceInput({inline: true, value: home, placeholder: 'Save to a space (optional)',
      tooltipText: SPACE_TOOLTIP});
    const mode = new ChoiceChip({value: 'move', items: [
      {value: 'move', label: 'Move', description: 'Move to the new space. The project gets its permissions; ' +
        'the previous location becomes a link.'},
      {value: 'link', label: 'Link', description: 'Reference the original project from the new space. ' +
        'You won\'t be able to edit it there without cloning first.'}]});
    const changed = computed(() => !same(input.value, home));
    const where = (t: DG.Project | null) => `"${SpaceInput.describe(t, input.tree.value)}"`;
    const hint = computed(() => {
      const target = input.value;
      if (home === null)
        return target === null ? '' : `Shared with whoever can see ${where(target)}.`;
      if (!changed.value)
        return `Lives in ${where(home)}.`;
      if (target === null)
        return `Moves out of ${where(home)} into your projects.`;
      return mode.signal.value === 'link' ? `Stays in ${where(home)}; ${where(target)} gets a link to it.` :
        `Moves from ${where(home)} to ${where(target)}.`;
    });
    self.effect(() => {
      const target = input.value;
      const text = projectName.value.trim();
      mode.root.hidden = !(home !== null && changed.value && target !== null);
      problem.value = null;
      // a link leaves the project where it is, so the platform checks no name for it
      if (target === null || !changed.value || text === '' || mode.signal.value === 'link')
        return;
      spaceHasName(target, text).then((taken) => {
        if (taken && input.peek() === target && projectName.peek().trim() === text) {
          problem.value = `${where(target)} already has a project named "${text}". ` +
            'Rename the project or pick another space.';
        }
      }).catch(() => {});
    });
    self.effect(() => {
      void input.value;
      void mode.signal.value;
      void problem.value;
      widget.onChanged();
    });
    return {
      root: div([divH([input, mode], {gap: 'm'}), statusLine(problem, hint)], {gap: 'xs'}),
      problem: () => problem.peek(),
      setName: (text: string) => projectName.value = text,
      target: () => input.peek(),
      setTarget: (t: DG.Project | null) => input.value = t,
      mode: () => mode.signal.peek() as 'move' | 'link',
      setMode: (m: 'move' | 'link') => mode.signal.value = m,
      resolveTarget: async (): Promise<SpaceSaveTarget> => {
        const target = input.peek();
        if (same(target, home))
          return {spaceId: home?.id ?? null, link: false, changed: false};
        if (target === null)
          return {spaceId: null, link: false, changed: true};
        rememberSpace(target);
        return {spaceId: target.id, link: home !== null && mode.signal.peek() === 'link', changed: true};
      },
    };
  }, {name: 'space-save-row'}));
  widget.root.append(row.root);
  Object.defineProperties(widget, {
    problem: {get: () => row.problem()},
    target: {get: () => row.target(), set: (t: DG.Project | null) => row.setTarget(t)},
    mode: {get: () => row.mode(), set: (m: 'move' | 'link') => row.setMode(m)},
  });
  widget.setName = (text) => row.setName(text);
  widget.resolveTarget = () => row.resolveTarget();
  return widget;
}

/** The row for a project owned by `homeId`: a space, or the user's own projects (none for a new project). */
export async function spaceSaveRowFor(homeId?: string | null, name?: string): Promise<SpaceSaveRow> {
  const home = homeId ? await grok.dapi.spaces.find(homeId).catch(() => null) : null;
  return spaceSaveRow(home, name ?? '');
}
