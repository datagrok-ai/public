import {BoolInput, Control, Form, NumberInput, batch, button, component, computed, effect, plural, signal, statusLine}
  from '@datagrok-libraries/u2';

export function orderForm(): Control {
  return component(() => {
    const quantity = signal<number | null>(10);
    const price = signal<number | null>(12.5);
    const member = signal(false);
    const total = computed(() => (quantity.value ?? 0) * (price.value ?? 0) * (member.value ? 0.9 : 1));
    const runs = signal(0);
    effect(() => {
      document.title = `Order: $${total.value.toFixed(2)}`;
      runs.value = runs.peek() + 1;
    });
    const form = new Form().addAll([
      new NumberInput({label: 'Quantity', mode: 'int', min: 1, bind: quantity}),
      new NumberInput({label: 'Unit price', min: 0, postfix: '$', bind: price}),
      new BoolInput({label: 'Member', switch: true, bind: member}),
    ]);
    const reset = () => batch(() => {
      quantity.value = 1;
      price.value = 10;
      member.value = false;
    });
    return form.addButtons((row) => row.append(button('Reset', reset),
      statusLine(form.shownProblem, computed(() => `Total $${total.value.toFixed(2)} · the effect ran ${plural(runs.value, 'time', 'times')}`))));
  }, {name: 'order-form'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = orderForm();
  host.append(demo.root);
  return demo;
}
