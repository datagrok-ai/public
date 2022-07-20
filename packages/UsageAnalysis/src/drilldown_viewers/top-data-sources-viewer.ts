import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {TopQueriesUsingDataSource} from "./top-queries-using-data-source";
import * as DG from "datagrok-api/dg";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {UaFilter} from "../filter2";
import {BehaviorSubject} from "rxjs"

export class TopDataSourcesViewer extends UaFilterableQueryViewer {
  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        'Data Sources',
        'TopDataSources',
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args) => {
            let viewer = new TopQueriesUsingDataSource(args.args.categories[0], filterStream);
            grok.shell.o = ui.block([viewer.root]);
          });
          if (t.rowCount > 0)
            return viewer.root
          else
            return ui.divText('Not enough data', {style:{color:'var(--red-3)', paddingBottom:'25px'}})
        }
    );
  }

}