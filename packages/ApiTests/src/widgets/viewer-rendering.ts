import {awaitCheck, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';

category('Viewer: rendering', () => {
  test('immediateRendering', async () => {
    const tv = grok.shell.addTableView(grok.data.demo.demog(100));
    const plot = tv.scatterPlot();
    let renders = 0;
    const sub = plot.onViewerRendered.subscribe(() => renders++);
    try {
      expect(plot.immediateRendering, false);
      await awaitCheck(() => renders > 0, 'plot did not paint', 3000);
      plot.immediateRendering = true;
      expect(plot.immediateRendering, true);
      renders = 0;
      plot.props.markerDefaultSize = 12;
      plot.props.markerOpacity = 50;
      await new Promise((r) => setTimeout(r, 0));
      expect(renders, 1);
    } finally {
      sub.unsubscribe();
      tv.close();
    }
  });
});
