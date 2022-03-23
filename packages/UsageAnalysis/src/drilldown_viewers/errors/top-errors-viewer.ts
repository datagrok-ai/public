import {UaFilterableQueryViewer} from "../../viewers/ua-filterable-query-viewer";
import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../../viewers/abstract/ua-query-viewer";
import {TopQueriesUsingDataSource} from "../top-queries-using-data-source";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import {PropertyPanel} from "../../property-panel";
import {UaDataFrameQueryViewer} from "../../viewers/ua-data-frame-query-viewer";
import {BehaviorSubject} from "rxjs"
import {ErrorMarkingPanel} from "../panels/error_marking_panel";

export class TopErrorsViewer extends UaFilterableQueryViewer {

  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        name,
        queryName,
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            let pp = new PropertyPanel(
                null,
                [new UaDataFrameQueryViewer(
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

            let errorMarkingPanel = new ErrorMarkingPanel();
            let accordion = await errorMarkingPanel.init(t, args.args.categories[0]);

            grok.shell.o = ui.divV([
              pp.getRoot(),
              accordion
            ]);

          });
          return viewer.root;
        }
    );
  }

}