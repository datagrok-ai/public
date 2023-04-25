import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilterableQueryViewer} from '../../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../../viewers/abstract/ua-query-viewer';
import {UaFilter} from '../../filter';
import {PropertyPanel} from '../../property-panel';
import {UaDataFrameQueryViewer} from '../../viewers/ua-data-frame-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class TopFunctionsViewer extends UaFilterableQueryViewer {
  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
      filterStream,
      'Functions',
      'TopFunctions',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args) => {
          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameQueryViewer(
              'Function Info',
              'FunctionInfoByName',
              (t: DG.DataFrame) => DG.Viewer.grid(t),
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            new UaDataFrameQueryViewer(
              'Users',
              'TopUsersOfFunction',
              (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions),
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            ],
            `Functions: ${args.args.categories[0]}`,
            'Functions');

          grok.shell.o = pp.getRoot();
        });
        return viewer;
      },
    );
  }
}
