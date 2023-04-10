import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilterableQueryViewer} from '../../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../../viewers/abstract/ua-query-viewer';
import {UaFilter} from '../../filter';
import {PropertyPanel} from '../../property-panel';
import {UaDataFrameQueryViewer} from '../../viewers/ua-data-frame-query-viewer';
import {BehaviorSubject} from 'rxjs';

export class TopErrorSourcesViewer extends UaFilterableQueryViewer {
  public constructor(filterStream: BehaviorSubject<UaFilter>) {
    super(
      filterStream,
      'Error Sources',
      'TopErrorSources',
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args) => {
          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameQueryViewer(
              'Function Info By Source',
              'FunctionInfoBySource',
              (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            ],
            `Error Sources: ${args.args.categories[0]}`,
            'Error Sources');

          grok.shell.o = pp.getRoot();
        });
        return viewer;
      },
    );
  }
}
