import {AsyncValue, ChoiceInput, Form, computed, divH, divV, span} from '@datagrok-libraries/u2';
import {readout} from './common';
import {sourceOf} from '../source-panel';
import {SRC_ROOT} from '../const';

export function chainedChoicesPage(): HTMLElement {
  return divV([
    span('Each choice loads its options from every pick above it, from a mock server that answers in 500 ms. ' +
      'Change a pick and everything below it reloads; a pick its new list lacks is cleared, which reloads what ' +
      'is below that. The code shown with the form is all of it.', 'u2demo-hint'),
    divH([chainedChoices(), divV([sourceOf(`${SRC_ROOT}/pages/chained-choices.ts`, 'chainedChoices')], 'u2demo-card')],
      'u2demo-funcs-ab'),
  ], 'u2demo-page u2demo-wide');
}

export function chainedChoices(): HTMLElement {
  // options load from every pick above, and again whenever one of them changes
  function dependentChoice(label: string, above: ChoiceInput[], marks: string) {
    const options = new AsyncValue(async () => {
      const picks = above.map((input) => input.peek());
      if (picks.includes(null))
        return [];
      await new Promise((resolve) => setTimeout(resolve, 500));
      return [...marks].map((mark) => `${picks[picks.length - 1]} ${mark}`);
    }, {deps: () => above.forEach((input) => input.value), debounceMs: 0});

    const choice = new ChoiceInput({label, items: [], nullable: true,
      loading: computed(() => options.state.value.kind === 'loading'),
      disabled: computed(() => {
        const empty = above.find((input) => input.value === null);
        return empty ? `Pick ${empty.label} first` : null;
      })});
    choice.effect(() => {
      const state = options.state.value;
      if (state.kind === 'ready')
        choice.setItems(state.value);
    });
    return choice;
  }

  const names = ['Ada', 'Ben', 'Cleo', 'Dan', 'Eve'];
  const chain = [new ChoiceInput({label: 'Name', items: names, nullable: true})];
  const marks = {Letter: 'abcde', Digit: '12345', Sign: '+-*/=', Shape: '○△□◇☆'};
  for (const [label, chars] of Object.entries(marks))
    chain.push(dependentChoice(label, [...chain], chars));

  const picked = computed(() => chain.map((input) => input.value ?? '–').join(' → '));
  return divV([new Form().addAll(chain).root, readout('picked', picked)]);
}
