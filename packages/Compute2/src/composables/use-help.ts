import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {Subject, from, of} from 'rxjs';
import {useSubscription} from '@vueuse/rxjs';
import {catchError, switchMap, tap} from 'rxjs/operators';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';

const helpInjectionKey = Symbol('HelpToken');

function useHelpService() {
  const helpHidden = Vue.ref(true);
  const helpLoading = Vue.ref(false);
  const helpContent = Vue.ref<string | undefined>(undefined);

  const func$ = new Subject<undefined | DG.Func>();

  const changeHelpFunc = (func: undefined | DG.Func) => {
    helpContent.value = undefined;
    func$.next(func ? Vue.markRaw(func) : undefined);
  };

  useSubscription(func$.pipe(
    tap(() => {
      helpLoading.value = true;
    }),
    switchMap((func: undefined | DG.Func) => {
      return (func && Utils.hasContextHelp(func)) ? from(Utils.getContextHelp(func)) : of(undefined);
    }),
    catchError((e) => {
      console.error(e);
      return of(undefined);
    }),
  ).subscribe((help) => {
    helpContent.value = help;
    helpLoading.value = false;
  }));

  return {
    helpHidden,
    helpContent,
    helpLoading,
    changeHelpFunc,
  };
}

export function setHelpService() {
  Vue.provide(helpInjectionKey, useHelpService());
}

export function useHelp() {
  const help = Vue.inject<ReturnType<typeof useHelpService> | undefined>(helpInjectionKey);
  return help ?? useHelpService();
}
