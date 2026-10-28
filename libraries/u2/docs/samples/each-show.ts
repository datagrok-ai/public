import {Control, NumberInput, button, component, computed, divH, divV, each, iconButton, show, signal, span}
  from '@datagrok-libraries/u2';
import './each-show.css';

interface Line {
  id: number;
  reagent: string;
  amount: number | null;
}

const LIMIT = 100;

export function reagentOrder(reagents: string[]): Control {
  return component(() => {
    let count = 0;
    const line = (): Line => {
      const id = ++count;
      return {id, reagent: reagents[(id - 1) % reagents.length], amount: 10};
    };
    const lines = signal<Line[]>([line(), line()]);
    const update = (id: number, amount: number | null) => lines.value = lines.value.map((l) => l.id === id ? {...l, amount} : l);
    const remove = (id: number) => lines.value = lines.value.filter((l) => l.id !== id);
    const total = computed(() => lines.value.reduce((sum, l) => sum + (l.amount ?? 0), 0));
    const rows = divV([], 'order-lines');
    each(rows, lines, (l) => String(l.id), (l) => divH([
      span(l.reagent, 'order-reagent'),
      new NumberInput({inline: true, name: `${l.reagent} amount`, value: l.amount, min: 0, postfix: 'mg',
        onChanged: (v) => update(l.id, v)}),
      iconButton('times', () => remove(l.id), {tooltip: `Remove ${l.reagent}`}),
    ], {gap: 'm'}), {update: () => true, empty: () => span('No reagents yet', 'order-note')});
    const warning = divH([]);
    show(warning, computed(() => total.value > LIMIT), () => span(`Over the ${LIMIT} mg a plate takes`, 'order-warning'));
    return [rows, divH([button('Add reagent', () => lines.value = [...lines.value, line()]),
      span(computed(() => `${lines.value.length} reagents, ${total.value} mg in all`), 'order-note'), warning], {gap: 'm'})];
  }, {name: 'reagent-order'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = reagentOrder(['DMSO', 'HEPES buffer', 'ATP', 'MgCl₂', 'DTT']);
  host.append(demo.root);
  return demo;
}
