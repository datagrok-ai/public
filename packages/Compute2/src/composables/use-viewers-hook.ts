import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useSubscription, from} from '@vueuse/rxjs'
import {BehaviorSubject, combineLatest, EMPTY, merge, Subject} from 'rxjs';
import {catchError, switchMap, tap, map, debounceTime, withLatestFrom, share, shareReplay} from 'rxjs/operators';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';

type ViewersData = Record<string, DG.Viewer | undefined>;

function runViewersHook(viewersHook: ViewersHook | undefined, io: string, viewers: ViewersData, meta: any) {
  if (viewersHook) {
    for (const [type, viewer] of Object.entries(viewers)) {
      if (viewer)
        viewersHook(io, type, viewer, meta);
    }
  }
}

export function useViewersHook(
  viewersHook: Vue.Ref<ViewersHook | undefined>,
  metaStates: Vue.Ref<Record<string, BehaviorSubject<any>> | undefined>,
  currentCall: Vue.Ref<DG.FuncCall>
) {
  const metaStates$ = from(metaStates, { immediate: true });
  const viewersHook$ = from(viewersHook, { immediate: true });
  const currentCall$ = from(currentCall, { immediate: true });

  const viewerStates$ = currentCall$.pipe(
    map((call) => {
      const nextStates: Record<string, BehaviorSubject<ViewersData>> = {};
      const states = [...Object.keys(call.inputs), ...Object.keys(call.outputs)];
      for (const io of states) {
        const viewers$ = new BehaviorSubject<ViewersData>({});
        nextStates[io] = viewers$;
      }
      return nextStates;
    }),
    shareReplay(1),
  );

  const viewerChanges$ = new Subject<readonly [DG.Viewer | undefined,  string, string]>();
  const viewersSub = viewerChanges$.pipe(
    withLatestFrom(viewerStates$),
    tap(([[viewer, io, type], viewerStates]) => {
      const viewers$ = viewerStates[io];
      if (viewers$) {
        const viewers = viewers$.value;
        viewers$.next({...viewers, [type]: viewer});
      }
    })
  ).subscribe();
  useSubscription(viewersSub);

  const viewHooksSub = combineLatest([viewersHook$, metaStates$, viewerStates$]).pipe(
    debounceTime(0),
    switchMap(([viewersHook, metaStates, viewerStates]) => {
      const hooksDeps: Record<string, { meta$?: BehaviorSubject<any>, viewers$?: BehaviorSubject<ViewersData> }> = {};
      for (const [io, meta$] of Object.entries(metaStates ?? {})) {
        if (!hooksDeps[io])
          hooksDeps[io] = {};
        hooksDeps[io].meta$ = meta$;
      }
      for (const [io, viewers$] of Object.entries(viewerStates ?? {})) {
        if (!hooksDeps[io])
          hooksDeps[io] = {};
        hooksDeps[io].viewers$ = viewers$;
      }
      const hookRuns = [...Object.entries(hooksDeps)].map(([io, {meta$, viewers$}]) => {
        if (!meta$ || !viewers$)
          return EMPTY;
        return combineLatest([viewers$, meta$]).pipe(
          tap(([viewers, meta]) => runViewersHook(viewersHook, io, viewers, meta)),
          catchError((e) => {
            grok.shell.error(e);
            console.error(e);
            return EMPTY;
          })
        )
      })
      return merge(...hookRuns);
    }),
  ).subscribe();
  useSubscription(viewHooksSub);

  const setViewerRef = (viewer: DG.Viewer | undefined, io: string, type: string) => {
    viewerChanges$.next([viewer, io, type] as const);
  }

  return { setViewerRef };
}
