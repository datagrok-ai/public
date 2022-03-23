import {UaFilterableQueryViewer} from "../../viewers/ua-filterable-query-viewer";
import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../../viewers/abstract/ua-query-viewer";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import {PropertyPanel} from "../../property-panel";
import {BehaviorSubject} from "rxjs"
import {UaDataFrameViewer} from "../../viewers/ua-data-frame-viewer";

export class TopFunctionErrorsViewer extends UaFilterableQueryViewer {

  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>, staticFilter?: Object, showName?: boolean) {
    super(
        filterStream,
        name,
        queryName,
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.setOptions({
            split: 'error_and_event',
          });
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            let df = t.rows.match({ error_and_event: args.args.categories[0]}).toDataFrame();

            let errorMessage = df.get('error_message', 0);
            let friendlyName = df.get('friendly_name', 0);

            let eventInfo = await grok.data.query('UsageAnalysis:EventByErrorMessageAndFriendlyName',
                {errorMessage: errorMessage, friendlyName: friendlyName});

            let pp = new PropertyPanel(
                null,
                [new UaDataFrameViewer(
                  'Errors Info',
                    eventInfo,
                  (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  false)],
            `Errors: ${args.args.categories[0]}`,
            'Errors');

            let isError = ui.boolInput('Is error', eventInfo.get('is_error', 0));
            let comment = ui.stringInput('Comment', eventInfo.get('comment', 0))

            let acc = DG.Accordion.create('ErrorInfo');
            acc.addPane('ErrorInfo', () => {

              let button = ui.buttonsInput([ui.bigButton('Save', async () => {
                await grok.data.query('UsageAnalysis:UpdateEventsIsErrorComment', {
                  errorMessage: errorMessage,
                  friendlyName: friendlyName,
                  isError: isError.value,
                  comment: comment.value
                });

                  grok.shell.info('Event type saved');
              })]);
              return ui.divV([isError.root, comment.root, button]);
            });


            grok.shell.o = ui.divV([
              pp.getRoot(),
              acc
            ]);

          });
          return viewer.root;
        },
        null,
        staticFilter,
        showName
    );
  }

}