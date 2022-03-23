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
import {TopFunctionErrorsViewer} from "../function_errors/top-function-errors-viewer";

export class TopPackagesViewer extends UaFilterableQueryViewer {

  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        name,
        queryName,
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            let entity = await grok.dapi.packages.filter(`shortName = "${args.args.categories[0]}"`).first();
            let pp = new PropertyPanel(
                entity,
                [new UaDataFrameQueryViewer(
                  'Package Info',
                  'PackageInfo',
                  (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              ),
              new UaDataFrameQueryViewer(
                  'Functions Of Package',
                  'TopFunctionsOfPackage',
                  (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              ),
              new UaDataFrameQueryViewer(
                  'Users Of Package',
                  'TopUsersOfPackage',
                  (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false
              ),
              new TopFunctionErrorsViewer('Errors Of Package','TopErrorsOfPackage', filterStream, {name: args.args.categories[0]}, false),
            ],
            `Packages: ${args.args.categories[0]}`,
            'Packages');

            grok.shell.o = pp.getRoot();
          });
          return viewer.root;
        }
    );
  }

}