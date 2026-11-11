import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import dayjs from 'dayjs';
import {getRunName, replaceForWindowsPath} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import type {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {showVueDialog} from '../../utils';
import {ExportOption} from '../../run-exports';
import {ExportRun, FailedExport, exportRuns, modelExports, modelKind, zipExports} from './bulk-export';
import {modelLabel, runCompleted} from './run-storage';
import './RunStorage.css';

type ModelExports = ExportOption[] | {error: string};

export const BulkExportDialog = Vue.defineComponent({
  name: 'BulkExportDialog',
  props: {
    locationName: {type: String, required: true},
    // the runs of the location and of everything below it, as the browser's filters keep them
    runs: {type: Array as Vue.PropType<ExportRun[]>, required: true},
    // the space whose own runs are exported without subspaces; none for My runs
    spaceId: {type: String, required: false},
    // the space and its subspaces that couldn't be listed
    unlisted: {type: Array as Vue.PropType<{spaceId: string, name: string}[]>, default: () => []},
    filtered: {type: Boolean, default: false},
    stepHistory: {type: Boolean, default: false},
    viewersHook: {type: Function as Vue.PropType<ViewersHook>, required: false},
  },
  emits: {
    done: (_result: null) => true,
  },
  setup(props, {emit}) {
    const includeSubspaces = Vue.ref(true);
    const exportsOf = Vue.shallowRef(new Map<string, ModelExports>());
    const chosen = Vue.ref<Record<string, string>>({});
    const progress = Vue.ref<{done: number, total: number}>();
    const failed = Vue.shallowRef<FailedExport[]>([]);
    const error = Vue.ref<string>();
    let isCancelled = false;

    const scoped = Vue.computed(() => includeSubspaces.value || !props.spaceId ? props.runs :
      props.runs.filter((item) => item.spaceId === props.spaceId));
    // runs saved as not completed are shown but not exported; the rest are checked once loaded
    const models = Vue.computed(() => {
      const groups = new Map<string, {model: DG.Func, items: ExportRun[], notCompleted: number}>();
      for (const item of scoped.value) {
        const group = groups.get(item.model.nqName) ?? {model: item.model, items: [], notCompleted: 0};
        if (runCompleted(item.run) === false)
          group.notCompleted++;
        else
          group.items.push(item);
        groups.set(item.model.nqName, group);
      }
      return [...groups.values()];
    });

    Vue.onMounted(() => {
      const all = new Map(props.runs.map((item) => [item.model.nqName, item.model]));
      for (const [nqName, model] of all) {
        modelExports(model, modelKind(model, props.stepHistory)).then(
          (options) => {
            exportsOf.value = new Map(exportsOf.value).set(nqName, options);
            if (options.length > 0)
              chosen.value = {...chosen.value, [nqName]: options[0].id};
          },
          (e) => exportsOf.value = new Map(exportsOf.value).set(nqName, {error: String(e?.message ?? e)}),
        );
      }
    });
    Vue.onUnmounted(() => isCancelled = true);

    const missing = Vue.computed(() => props.unlisted.filter(({spaceId}) =>
      includeSubspaces.value || spaceId === props.spaceId).map(({name}) => name));
    const isLoading = Vue.computed(() => models.value.some(({model}) => !exportsOf.value.has(model.nqName)));
    const jobs = Vue.computed(() => models.value.flatMap(({model, items}) => {
      const exportId = chosen.value[model.nqName];
      const kind = modelKind(model, props.stepHistory);
      return exportId ? items.map((item) => ({item, kind, exportId})) : [];
    }));

    const start = async () => {
      if (isLoading.value || jobs.value.length === 0 || progress.value)
        return;
      progress.value = {done: 0, total: jobs.value.length};
      try {
        const result = await exportRuns(jobs.value, {
          onProgress: (done) => progress.value = {done, total: progress.value!.total},
          isCancelled: () => isCancelled,
          viewersHook: props.viewersHook,
        });
        if (isCancelled)
          return;
        if (result.files.length > 0) {
          const name = `${replaceForWindowsPath(props.locationName)} - ${dayjs().format('YYYY-MM-DD')}.zip`;
          DG.Utils.download(name, new Blob([zipExports(result.files) as any]));
        }
        if (result.notCompleted.length > 0) {
          grok.shell.warning(`Skipped ${result.notCompleted.length} runs that aren't completed: ` +
            result.notCompleted.map(getRunName).join(', '));
        }
        if (result.failed.length === 0)
          emit('done', null);
        failed.value = result.failed;
      } catch (e: any) {
        error.value = String(e?.message ?? e);
        progress.value = undefined;
      }
    };
    const cancel = () => {
      isCancelled = true;
      emit('done', null);
    };

    const renderChoice = (model: DG.Func) => {
      const options = exportsOf.value.get(model.nqName);
      if (!options)
        return <span class='c2-runs-muted'>Loading…</span>;
      if ('error' in options)
        return <span class='c2-runs-error' title={options.error}>Couldn't load the exports</span>;
      if (options.length === 0)
        return <span class='c2-runs-muted'>No exports</span>;
      if (options.length === 1)
        return <span>{options[0].name}</span>;
      return <select class='ui-input-editor' value={chosen.value[model.nqName]} disabled={!!progress.value}
        onChange={(e: Event) =>
          chosen.value = {...chosen.value, [model.nqName]: (e.target as HTMLSelectElement).value}}>
        { options.map(({id, name}) => <option key={id} value={id}>{name}</option>) }
      </select>;
    };

    const renderFailed = () => <div class='c2-export-failed'>
      <div class='c2-runs-error'>{`Couldn't export ${failed.value.length} of ${progress.value!.total} runs:`}</div>
      { failed.value.map(({run, error}) => <div key={run.id}>
        <span class='c2-runs-row-name'>{getRunName(run)}</span>: <span class='c2-runs-muted'>{error}</span>
      </div>) }
    </div>;

    const footerNote = () => {
      if (progress.value)
        return failed.value.length > 0 ? '' : `Exporting ${progress.value.done} of ${progress.value.total} runs…`;
      return `${jobs.value.length} ${jobs.value.length === 1 ? 'run' : 'runs'} to export`;
    };

    return () => (
      <div class='c2-runs-dialog c2-export-dialog'>
        { props.spaceId &&
          <label class='c2-export-toggle'>
            <input type='checkbox' checked={includeSubspaces.value} disabled={!!progress.value}
              onChange={(e: Event) => includeSubspaces.value = (e.target as HTMLInputElement).checked}/>
            Include subspaces
          </label> }
        { missing.value.length > 0 &&
          <div class='c2-export-warning'>
            <i class='fas fa-exclamation-triangle'/>
            {`Couldn't list ${missing.value.join(', ')}, so ${missing.value.length === 1 ? 'its' : 'their'} ` +
              'runs aren\'t exported.'}
          </div> }
        { props.filtered &&
          <div class='c2-export-warning'>
            <i class='fas fa-exclamation-triangle'/>
            Search or filters are on: only the runs they show are exported.
          </div> }
        <div class='c2-export-models'>
          { models.value.map(({model, items, notCompleted}) => {
            const count = items.length + notCompleted;
            const dimmed = items.length === 0 ? 'c2-runs-dimmed' : '';
            return [
              <span key={`${model.nqName}-name`} class={['c2-runs-row-name', dimmed]}>{modelLabel(model)}</span>,
              <span key={`${model.nqName}-count`} class={['c2-runs-muted', dimmed]}>
                {`${count} ${count === 1 ? 'run' : 'runs'}`}
                {notCompleted > 0 && items.length > 0 && `, ${notCompleted} not completed`}</span>,
              <span key={`${model.nqName}-export`} class={dimmed}>
                {items.length === 0 ? 'Not completed' : renderChoice(model)}</span>,
            ];
          }) }
          { models.value.length === 0 && <span class='c2-runs-muted'>No runs to export</span> }
        </div>
        { failed.value.length > 0 && renderFailed() }
        { error.value && <div class='c2-runs-error'>{error.value}</div> }
        <div class='c2-runs-buttons'>
          <span class='c2-runs-muted c2-runs-footer-note'>{footerNote()}</span>
          <div class='c2-runs-button-group'>
            { failed.value.length > 0 ?
              <button class='ui-btn ui-btn-ok' onClick={() => emit('done', null)}>Close</button> :
              <>
                <button class='ui-btn' onClick={cancel}>Cancel</button>
                <button class={['ui-btn', 'ui-btn-ok',
                  isLoading.value || jobs.value.length === 0 || progress.value ? 'd4-disabled' : '']}
                onClick={start}>Export</button>
              </> }
          </div>
        </div>
      </div>
    );
  },
});

/** `modelsName` is the browser's model filter: a model or All models. */
export function openBulkExportDialog({modelsName, ...props}: {
  locationName: string, modelsName: string, runs: ExportRun[], spaceId?: string,
  unlisted: {spaceId: string, name: string}[], filtered: boolean, stepHistory: boolean, viewersHook?: ViewersHook,
}) {
  return showVueDialog<null>(`Export runs of ${props.locationName} · ${modelsName}`, BulkExportDialog, props,
    {width: 620, height: 420});
}
