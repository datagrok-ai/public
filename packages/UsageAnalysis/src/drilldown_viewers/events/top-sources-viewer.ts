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

export class TopSourcesViewer extends UaFilterableQueryViewer {

  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        'Sources',
        'TopSources',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args) => {

            let pp = new PropertyPanel(
                null,
                null,
                [new UaDataFrameQueryViewer(
                  'Functions',
                  'TopFunctionsOfSource',
                  (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              ),
              new UaDataFrameQueryViewer(
                  'Users',
                  'TopUsersOfSource',
                  (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              )
            ],
            `Sources: ${args.args.categories[0]}`,
            'Sources');

            grok.shell.o = pp.getRoot();
          });
          if (t.rowCount > 0)
            return viewer.root
          else
            return ui.divText('Not enough data', {style:{color:'var(--red-3)', paddingBottom:'25px'}})
        }
    );
  }

}