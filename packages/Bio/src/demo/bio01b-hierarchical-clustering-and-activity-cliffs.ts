import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {awaitCheck} from '@datagrok-libraries/test/src/test';

export async function demoActivityCliffsCyclicLayout(): Promise<void> {
  grok.shell.windows.showContextPanel = true;
  const p = await grok.functions.eval('Bio:BioDemoActivityCliffs');
  const project = await grok.dapi.projects.find(p.id);
  await project.open();
  let scatterPlot: DG.Viewer | null = null;
  for (const i of grok.shell.tv.viewers) {
    if (i.type == DG.VIEWER.SCATTER_PLOT)
      scatterPlot = i;
  }
  let cliffsLink;
  try {
    await awaitCheck(() => {
      const link = scatterPlot?.root.getElementsByClassName('scatter_plot_link');
      if (link?.length) {
        cliffsLink = link[0];
        return true;
      }
      return false;
    }, '', 10000);
    (cliffsLink as any as HTMLElement).click();
  } catch (e) {}
}
