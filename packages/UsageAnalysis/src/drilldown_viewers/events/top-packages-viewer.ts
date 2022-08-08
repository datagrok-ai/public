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
import { ViewHandler } from "../../view-handler";
import {TopPackageFunctionsViewer} from "./top-package-functions-viewer";

export class TopPackagesViewer extends UaFilterableQueryViewer {

  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
        filterStream,
        name,
        queryName,
        (t: DG.DataFrame) => {
          let viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
          viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
            this.categorySelected(args.args.categories[0]);
          });
          return viewer.root;
        }
    );    
  }

  async categorySelected(category: string)  {
    ViewHandler.getInstance().setUrlParam('package', category);

    let entity = await grok.dapi.packages.filter(`shortName = "${category}"`).first();
    let pp = new PropertyPanel(
        entity,
        new UaDataFrameQueryViewer(
            'Package Info',
            'PackageInfo',
            (t: DG.DataFrame) => {
              let res: any = {};
              for (let c of t.columns) {
                res[c.name] = c.get(0);
              }

              return ui.tableFromMap(res);
            },
            null as any,
            {name: category},
            this.filterSubscription.getValue(),
            false
        ),
        [
      new TopPackageFunctionsViewer(
          'Functions',
          'TopFunctionsOfPackage',
          this.filterSubscription,
          category,
          false
      ),
      new UaDataFrameQueryViewer(
          'Users',
          'TopUsersOfPackage',
          (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions).root,
          null as any,
          {name: category},
          this.filterSubscription.getValue(),
          false
      ),
      new TopFunctionErrorsViewer('Errors','TopErrorsOfPackage', this.filterSubscription, {name: category}, false),
    ],
    `Packages: ${category}`,
    'Packages');

    grok.shell.o = pp.getRoot();
  }

}