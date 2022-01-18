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

export class TopErrorsViewer extends UaFilterableViewer {

  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        'Errors',
        'TopErrors',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            let pp = new PropertyPanel(
                null,
                [new UaDataFrameViewer(
                  'Errors Info',
                  'FunctionInfoByFriendlyName',
                  (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              )
            ],
            `Errors: ${args.args.categories[0]}`,
            'Errors');

            let id = t.rows.match({ friendly_name: args.args.categories[0]}).toDataFrame().get('id', 0);

            let error = await grok.dapi.logTypes.include('message,isError').filter(`id = "${id}"`).first();

            let isError = ui.boolInput('Is error', error.isError);
            let comment = ui.stringInput('Comment', error.comment)

            let acc = DG.Accordion.create('ErrorInfo');
            acc.addPane('ErrorInfo', () => {

              let button = ui.buttonsInput([ui.bigButton('Save', () => {
                error.isError = isError.value;
                error.comment = comment.value;
                grok.dapi.logTypes.save(error);
                grok.shell.info('Event type saved');
              })])
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