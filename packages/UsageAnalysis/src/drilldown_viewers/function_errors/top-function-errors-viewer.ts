import {UaFilterableViewer} from "../../viewers/ua-filterable-viewer";
import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../../viewers/ua-query-viewer";
import {TopQueriesUsingDataSource} from "../top-queries-using-data-source";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import {PropertyPanel} from "../../property-panel";
import {UaDataFrameViewer} from "../../viewers/ua-data-frame-viewer";
import {BehaviorSubject} from "rxjs"
import {ErrorMarkingPanel} from "../panels/error_marking_panel";

export class TopFunctionErrorsViewer extends UaFilterableViewer {

  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>) {
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

            let pp = new PropertyPanel(
                null,
                [new UaDataFrameViewer(
                  'Errors Info',
                  'EventByErrorMessageAndFriendlyName',
                  (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  null as any,
                  {errorMessage: errorMessage, friendlyName: friendlyName},
                  filterStream.getValue(),
                  false
              )
            ],
            `Errors: ${args.args.categories[0]}`,
            'Errors');

            let isError = ui.boolInput('Is error', df.get('is_error', 0));
            let comment = ui.stringInput('Comment', df.get('comment', 0))

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
        }
    );
  }

}