import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilterableQueryViewer} from '../../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../../viewers/abstract/ua-query-viewer';
import {UaFilter} from '../../filter';
import {PropertyPanel} from '../../property-panel';
import {UaDataFrameQueryViewer} from '../../viewers/ua-data-frame-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class TopSourcesViewer extends UaFilterableQueryViewer {
  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
      filterStream,
      'Sources',
      'TopSources',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args) => {
          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameQueryViewer(
              'Functions',
              'TopFunctionsOfSource',
              (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions),
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            new UaDataFrameQueryViewer(
              'Users',
              'TopUsersOfSource',
              (t: DG.DataFrame) => DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions),
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            ],
            `Sources: ${args.args.categories[0]}`,
            'Sources');

          grok.shell.o = pp.getRoot();
        });
        return viewer;
      },
    );
  }
}
