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

export class TopUsersViewer extends UaFilterableQueryViewer {

  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        'Users',
        'TopUsers',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            let entity = await grok.dapi.users.filter(`login = "${args.args.categories[0]}"`).first();
            let pp = new PropertyPanel(
                entity,
                [new UaDataFrameQueryViewer(
                  'Events Of User',
                  'TopEventsOfUser',
                  (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              ),

            ],
            `Users: ${args.args.categories[0]}`,
            'Users');

            grok.shell.o = pp.getRoot();
          });
          return viewer.root;
        }
    );
  }

}