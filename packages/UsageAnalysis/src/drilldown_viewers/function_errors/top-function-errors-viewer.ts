/*
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilterableQueryViewer} from '../../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../../viewers/abstract/ua-query-viewer';
import {UaFilter} from '../../filter';
import {PropertyPanel} from '../../property-panel';
import {BehaviorSubject} from 'rxjs';
import {UaDataFrameViewer} from '../../viewers/ua-data-frame-viewer';

export class TopFunctionErrorsViewer extends UaFilterableQueryViewer {
  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>, staticFilter?: Object,
    showName?: boolean) {
    super(
      filterStream,
      name,
      queryName,
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
        viewer.setOptions({
          split: 'error_and_event',
        });
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          const df = t.rows.match({error_and_event: args.args.categories[0]}).toDataFrame();

          const errorMessage = df.get('error_message', 0);
          const friendlyName = df.get('friendly_name', 0);

          const eventInfo = await grok.data.query('UsageAnalysis:EventByErrorMessageAndFriendlyName',
            {errorMessage: errorMessage, friendlyName: friendlyName});

          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameViewer(
              'Errors',
              eventInfo,
              (t: DG.DataFrame) => {
                const grid = DG.Viewer.grid(t);
                return grid.root;
              },
            )],
            `Errors: ${args.args.categories[0]}`,
            'Errors');

          const isError = ui.input.bool('Is error', {value: eventInfo.get('is_error', 0)});
            const comment = ui.input.string('Comment', {value: eventInfo.get('comment', 0)});
          const acc = DG.Accordion.create('ErrorInfo');

          acc.addPane('Error details', () => {
            const button = ui.buttonsInput([ui.bigButton('Save', async () => {
              await grok.data.query('UsageAnalysis:UpdateEventsIsErrorComment', {
                errorMessage: errorMessage,
                friendlyName: friendlyName,
                isError: isError.value,
                comment: comment.value,
              });
              grok.shell.info('Event type saved');
            })]);
            return ui.divV([ui.h2(eventInfo.get('error_message', 0)), ui.divText(eventInfo.get('error_stack_trace', 0),
              {style: {maxHeight: '250px', overflowY: 'scroll'}}), isError.root, comment.root, button]);
          });
          grok.shell.o = ui.divV([
            pp.getRoot(),
            acc,
          ]);
        });
        return viewer;
      },
      null,
      staticFilter,
    );
  }
}
*/
