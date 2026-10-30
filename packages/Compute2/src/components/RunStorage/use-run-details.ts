import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {EMPTY, Subject, from, of} from 'rxjs';
import {catchError, switchMap} from 'rxjs/operators';
import {useSubscription} from '@vueuse/rxjs';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {getRunAuthorId} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {Place, canEditSpace} from './run-storage';

export type RunDetails = {run: DG.FuncCall, author: string, modifiedBy?: string, canEdit: boolean, canUnlink: boolean};

/** The run shown in the details pane: loaded with its authors, and whether the user may edit it or take it out of
 * the space it is in. */
export function useRunDetails(showDetails: () => boolean, spaceOf: (run: DG.FuncCall) => Place) {
  const current = Vue.shallowRef<DG.FuncCall>();
  const details = Vue.shallowRef<RunDetails>();
  const isLoadingDetails = Vue.ref(false);
  const mode = Vue.ref<'view' | 'edit'>('view');
  const actionError = Vue.ref<string>();

  const loadDetails = async (run: DG.FuncCall): Promise<RunDetails> => {
    const [loaded, editSpace] = await Promise.all([
      historyUtils.shallowLoadRun(run.id),
      canEditSpace(spaceOf(run)),
    ]);
    // after someone else edited the details, the platform author is the last editor
    const authorId = getRunAuthorId(loaded);
    const lastEditor = loaded.author;
    const author = authorId && authorId !== lastEditor?.id ?
      await grok.dapi.users.find(authorId).then((user) => user?.friendlyName ?? '', () => '') : undefined;
    return {
      run: Vue.markRaw(loaded),
      author: author ?? lastEditor?.friendlyName ?? '',
      modifiedBy: author !== undefined ? lastEditor?.friendlyName : undefined,
      canEdit: historyUtils.canEditRunDetails(authorId, editSpace),
      canUnlink: editSpace,
    };
  };

  const chosen$ = new Subject<DG.FuncCall | undefined>();
  useSubscription(chosen$.pipe(
    switchMap((run) => {
      isLoadingDetails.value = !!run && showDetails();
      if (!run || !isLoadingDetails.value)
        return EMPTY;
      return from(loadDetails(run)).pipe(catchError((e) => {
        actionError.value = `Couldn't load the run: ${e?.message ?? e}`;
        return of(undefined);
      }));
    }),
  ).subscribe((loaded) => {
    details.value = loaded;
    isLoadingDetails.value = false;
  }));

  const choose = (run: DG.FuncCall | undefined) => {
    current.value = run;
    details.value = undefined;
    mode.value = 'view';
    actionError.value = undefined;
    chosen$.next(run);
  };

  const startEdit = () => {
    actionError.value = undefined;
    mode.value = 'edit';
  };
  const stopEdit = () => mode.value = 'view';

  return {
    current: current as Readonly<typeof current>,
    details: details as Readonly<typeof details>,
    isLoadingDetails: isLoadingDetails as Readonly<typeof isLoadingDetails>,
    mode: mode as Readonly<typeof mode>,
    actionError: actionError as Readonly<typeof actionError>,
    choose, startEdit, stopEdit,
  };
}
