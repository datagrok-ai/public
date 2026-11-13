import * as Vue from 'vue';
import {BehaviorSubject} from 'rxjs';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {useUnwrappedCallMeta} from '@datagrok-libraries/webcomponents-vue/src/composables/useUnwrappedCallMeta';

function setup(initial?: Record<string, BehaviorSubject<any>>) {
  const source = Vue.shallowRef(initial);
  const scope = Vue.effectScope();
  const state = scope.run(() => useUnwrappedCallMeta(() => source.value))!;
  return {source, state, stop: () => scope.stop()};
}

category('WebComponentsVue: useUnwrappedCallMeta', () => {
  test('unwraps the meta subjects and follows them', async () => {
    const a = new BehaviorSubject<any>({hidden: true});
    const {state, stop} = setup({a});
    try {
      expectDeepEqual({...state}, {a: {hidden: true}});
      a.next({items: ['x']});
      expectDeepEqual({...state}, {a: {items: ['x']}});
    } finally {
      stop();
    }
  });

  test('a new source drops the keys of the old one', async () => {
    const a = new BehaviorSubject<any>({hidden: true});
    const {source, state, stop} = setup({a});
    try {
      source.value = {b: new BehaviorSubject<any>({hidden: false})};
      await Vue.nextTick();
      expectDeepEqual({...state}, {b: {hidden: false}});
      a.next({hidden: false});
      expectDeepEqual({...state}, {b: {hidden: false}}, {prefix: 'the old subject still writes'});
    } finally {
      stop();
    }
  });

  test('a cleared source drops the old keys', async () => {
    const {source, state, stop} = setup({a: new BehaviorSubject<any>({hidden: true})});
    try {
      source.value = undefined;
      await Vue.nextTick();
      expectDeepEqual({...state}, {}, {prefix: 'the meta of the previous call stayed'});
    } finally {
      stop();
    }
  });
});
