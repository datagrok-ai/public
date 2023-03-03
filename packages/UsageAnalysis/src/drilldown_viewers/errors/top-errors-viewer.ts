import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {UaFilterableQueryViewer} from '../../viewers/ua-filterable-query-viewer';
import {UaQueryViewer} from '../../viewers/abstract/ua-query-viewer';
import {UaFilter} from '../../filter';
import {PropertyPanel} from '../../property-panel';
import {UaDataFrameQueryViewer} from '../../viewers/ua-data-frame-query-viewer';
import {BehaviorSubject} from 'rxjs';
import {ErrorMarkingPanel} from '../panels/error_marking_panel';

export class TopErrorsViewer extends UaFilterableQueryViewer {
  public constructor(name: string, queryName: string, filterStream: BehaviorSubject<UaFilter>) {
    super(
      filterStream,
      name,
      queryName,
      (t: DG.DataFrame) => {
        const viewer = DG.Viewer.barChart(t, UaQueryViewer.defaultBarchartOptions);
        viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe(async (args) => {
          const pp = new PropertyPanel(
            null,
            null,
            [new UaDataFrameQueryViewer(
              'Errors Info',
              'FunctionInfoByFriendlyName',
              (t: DG.DataFrame) => DG.Viewer.grid(t).root,
                  null as any,
                  {name: args.args.categories[0]},
                  filterStream.getValue(),
                  false,
            ),
            ],
            `Errors: ${args.args.categories[0]}`,
            'Errors');

          const errorMarkingPanel = new ErrorMarkingPanel();
          const accordion = await errorMarkingPanel.init(t, args.args.categories[0]);

          grok.shell.o = ui.divV([
            pp.getRoot(),
            accordion,
          ]);
        });
        return viewer.root;
      },
    );
  }
}
