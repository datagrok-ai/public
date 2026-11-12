import {Control, TextInput, component, computed, h2} from '@datagrok-libraries/u2';

export interface Greeting {
  readonly person: TextInput;
}

export function greeting(): Control & Greeting {
  return component(() => {
    const person = new TextInput({label: 'Your name', value: 'Ada'});
    const text = computed(() => `Hello, ${person.value.trim() || 'stranger'}!`);
    return {root: [h2(text), person], person};
  }, {name: 'greeting'});
}

// Live demo
export function render(host: HTMLElement): Control {
  const demo = greeting();
  host.append(demo.root);
  return demo;
}
