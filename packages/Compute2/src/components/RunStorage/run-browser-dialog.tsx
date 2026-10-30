import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {RunBrowser} from './RunBrowser';
import {DialogSize, showVueDialog} from '../../utils';
import {isRoleOnlyModelFilter} from '../../package-instance';

const RunBrowserDialog = Vue.defineComponent({
  name: 'RunBrowserDialog',
  props: {
    model: {type: DG.Func, required: false},
    selectModel: {type: Boolean, default: false},
    allowAllModels: {type: Boolean, default: false},
    stepHistory: {type: Boolean, default: false},
    openSeveral: {type: Boolean, default: false},
  },
  emits: {done: (_runs: DG.FuncCall[] | null) => true},
  setup(props, {emit}) {
    return () => <RunBrowser model={props.model} selectModel={props.selectModel}
      allowAllModels={props.allowAllModels} stepHistory={props.stepHistory} openSeveral={props.openSeveral}
      roleOnlyFilter={isRoleOnlyModelFilter()} layout='dialog'
      onOpen={(run: DG.FuncCall) => emit('done', [run])}
      onOpenSeveral={(runs: DG.FuncCall[]) => emit('done', runs)} onCancel={() => emit('done', null)}/>;
  },
});

/** The run browser as a dialog; resolves the runs to open (one unless openSeveral), or null. */
export function openRunBrowserDialog(options: {
  title: string, model?: DG.Func, selectModel?: boolean, allowAllModels?: boolean, stepHistory?: boolean,
  openSeveral?: boolean, size?: DialogSize,
}) {
  const {title, size, ...props} = options;
  return showVueDialog<DG.FuncCall[]>(title, RunBrowserDialog, props, size ?? {width: 980, height: 600});
}

/** Model Hub's Browse Runs: runs of every model, or of the one chosen, opened in the model's own view. */
export async function browseRuns() {
  // a blocking dialog of its own, so it takes most of the window
  const size = {width: Math.round(window.innerWidth * 0.9), height: Math.round(window.innerHeight * 0.9),
    fullScreen: true};
  const runs = await openRunBrowserDialog({title: 'Browse runs', selectModel: true, allowAllModels: true,
    openSeveral: true, size});
  for (const run of runs ?? []) {
    try {
      await grok.functions.call('Compute2:OpenWorkflowRun', {id: run.id});
    } catch (e: any) {
      grok.shell.error(`Couldn't open ${run.friendlyName}: ${e?.message ?? e}`);
    }
  }
}
