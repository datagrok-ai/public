import {ChoiceChip, Control, component, computed, divH, each, signal, statusLine} from '@datagrok-libraries/u2';
import type {ChoiceChipItem} from '@datagrok-libraries/u2';

const CAPS: ChoiceChipItem[] = [
  {value: 'H', label: 'H', description: 'hydrogen'},
  {value: 'OH', label: 'OH', description: 'hydroxyl'},
  {value: 'Br', label: 'Br', description: 'bromo'},
  {value: 'X', label: 'X', description: 'a placeholder, not a cap'},
];

export function capsEditor(labels: string[]): Control {
  return component(() => {
    const caps = signal(new Map(labels.map((label) => [label, label === 'R1' ? 'H' : 'OH'])));
    const problem = (cap: string | undefined) => cap === 'X' ? 'X is a placeholder: pick the cap it gets' : null;
    const chips = divH([], {gap: 's'});
    each(chips, () => labels, (label) => label, (label) => new ChoiceChip({label, items: CAPS,
      value: caps.peek().get(label), invalid: computed(() => problem(caps.value.get(label))),
      onInput: (cap) => caps.value = new Map(caps.peek()).set(label, cap ?? '')}).root);
    const first = computed(() => labels.map((l) => problem(caps.value.get(l))).find((p) => p !== null) ?? null);
    return [chips, statusLine(first, computed(() => labels.map((l) => `${l} ${caps.value.get(l)}`).join(', ')))];
  }, {name: 'caps-editor'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const editor = capsEditor(['R1', 'R2', 'R3']);
  host.append(editor.root);
  return editor;
}
