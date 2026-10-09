import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {RibbonMenuItem} from '@datagrok-libraries/webcomponents-vue';

/** Report a bug / Request a feature items from the REPORT_BUG_URL and REQUEST_FEATURE_URL package settings. */
export function useFeedbackItems(func: Vue.Ref<DG.Func | undefined>) {
  return Vue.computed<RibbonMenuItem[]>(() => {
    const settings = func.value?.package?.settings;
    const items: RibbonMenuItem[] = [];
    if (settings?.REPORT_BUG_URL)
      items.push({text: 'Report a bug', onClick: () => window.open(settings.REPORT_BUG_URL, '_blank')});
    if (settings?.REQUEST_FEATURE_URL)
      items.push({text: 'Request a feature', onClick: () => window.open(settings.REQUEST_FEATURE_URL, '_blank')});
    return items;
  });
}
