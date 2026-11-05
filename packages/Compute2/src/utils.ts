import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {dfToViewerMapping} from '@datagrok-libraries/compute-utils/shared-utils/export';
import {DEFAULT_FLOAT_FORMAT} from '@datagrok-libraries/webcomponents-vue';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {BehaviorSubject} from 'rxjs';
import * as Vue from 'vue';

declare let ENABLE_VUE_DEV_TOOLS: any;

export function setVueAppOptions(app: Vue.App<any>) {
  app.config.compilerOptions.isCustomElement = (tag) => tag.startsWith('dg-') || tag === 'dock-spawn-ts';
  if (ENABLE_VUE_DEV_TOOLS)
    app.config.performance = true;
}

// fullScreen lifts the platform's cap on dialog content height, so a large dialog fills its size
export type DialogSize = {width: number, height: number, fullScreen?: boolean};

/** Mounts a Vue component into a modal platform dialog and resolves what the component reports. */
export function showVueDialog<T>(
  title: string,
  component: Vue.Component,
  props: Record<string, any>,
  size: DialogSize,
): Promise<T | null> {
  return new Promise((resolve) => {
    const host = ui.div([], {style: {width: '100%', height: '100%'}});
    // the component renders its own buttons
    const dialog = ui.dialog({title, showFooter: false});
    let settled = false;
    const app = Vue.createApp(component, {
      ...props,
      onDone: (result: T | null) => {
        settled = true;
        resolve(result);
        dialog.close();
      },
    });
    setVueAppOptions(app);
    dialog.add(host);
    dialog.onClose.subscribe(() => {
      app.unmount();
      if (!settled)
        resolve(null);
    });
    dialog.show({modal: true, resizable: true, ...size, center: true});
    app.mount(host);
  });
}

// export-time viewers are created ad hoc and never mounted; detach releases their dart side
export function disposeViewers(mapping: {[key: string]: (DG.Viewer | undefined)[]} | undefined) {
  for (const viewers of Object.values(mapping ?? {})) {
    for (const viewer of viewers)
      viewer?.detach();
  }
}

export async function getViewers(call: DG.FuncCall, viewersHook?: ViewersHook, metaState?: Record<string, BehaviorSubject<any>>) {
  const mappings = await dfToViewerMapping(call);
  if (viewersHook) {
    for (const [ioName, viewers] of Object.entries(mappings ?? {})) {
      for (const viewer of viewers) {
        if (!viewer)
          continue;
        const meta = metaState?.[ioName]?.value;
        viewersHook(ioName, viewer.type, viewer, meta);
      }
    }
  }
  return mappings;
}

export function applyDefaultGridFloatFormat(viewer: DG.Viewer | undefined, type: string) {
  if (!viewer || type !== DG.VIEWER.GRID) return;
  const grid = viewer as DG.Grid;
  for (let i = 0; i < grid.columns.length; i++) {
    const gc = grid.columns.byIndex(i);
    const col = gc?.column;
    if (!col || col.type !== DG.COLUMN_TYPE.FLOAT) continue;
    if (gc!.format) continue;
    if (col.tags?.['format']) continue;
    gc!.format = DEFAULT_FLOAT_FORMAT;
  }
}

// Guards result-dependent actions (export/save) shared by RFV and RFVApp. Returns true when the
// action may proceed; otherwise surfaces a shell message explaining why the results aren't ready.
export function canUseResults(
  state: {isRunning?: boolean, runError?: unknown, isOutputOutdated?: boolean} | undefined,
  action: string,
): boolean {
  if (state?.isRunning) {
    grok.shell.warning(`The model is still running — wait for it to finish before ${action}.`);
    return false;
  }
  if (state?.runError) {
    grok.shell.error(`The last run finished with an error — fix the inputs and rerun before ${action}.`);
    return false;
  }
  if (state?.isOutputOutdated) {
    grok.shell.warning(`Results are outdated — run the model to update them before ${action}.`);
    return false;
  }
  return true;
}

export function pinView(view?: DG.ViewBase): void {
  if (view && !view.isPinned)
    view.pin();
}

// shared inline-style colors; tailwind arbitrary-value classes stay literal (the JIT needs them static)
export const STICKY_BAR_BACKGROUND = 'color-mix(in srgb, var(--white) 75%, transparent)';
