import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';
import {loadUsers, setupUserIconRenderer} from '../utils';
import '../../css/usage_analysis.css';

const filters = ui.box();
filters.classList.add('ua-filters');

const filtersStyle = {
  columnNames: ['source', 'user', 'event_time'],
};

const sourceColors: {[key: string]: string} = {
  'error': '#d62727',
  'audit': '#1f77b4',
  'data-query': '#9467bd',
  'function': '#2ba02b',
  'function-package': '#97df8a',
  'script': '#ffbb78',
};

export class LogView extends UaView {
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Log';
  }

  async initViewers(path?: string): Promise<void> {
    const users = await loadUsers();

    const logViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Log',
      queryName: 'LogTail',
      createViewer: (t: DG.DataFrame) => {
        const viewer = DG.Viewer.grid(t, {
          'showColumnLabels': false,
          'showRowHeader': false,
          'showColumnGridlines': false,
          'allowRowSelection': false,
          'allowBlockSelection': false,
          'showCurrentCellOutline': false,
        });
        filters.append(DG.Viewer.filters(t, filtersStyle).root);

        viewer.columns.setOrder(['source', 'user', 'event_time', 'event_description', 'id']);
        viewer.col('source')!.width = 30;
        viewer.col('ugid')!.visible = false;

        viewer.onCellPrepare(function(gc) {
          if (gc.gridColumn.name === 'event_description')
            gc.style.font = '13px monospace';

          if (gc.gridColumn.name === 'event_time' || gc.gridColumn.name === 'id') {
            gc.style.textColor = 0xFFB8BAC0;
            gc.style.font = '13px Roboto';
          }
        });

        setupUserIconRenderer(viewer, users, ['user']);

        viewer.onCellRender.subscribe(function(args) {
          if (args.cell.gridColumn.name == 'source') {
            args.g.beginPath();
            args.g.fillStyle = sourceColors[args.cell.cell.value] ?? '#7f7f7f';
            args.g.fillRect(args.bounds.x + 12, args.bounds.y + 11, 8, 8);
            args.preventDefault();
          }
        });
        return viewer;
      },
    });

    logViewer.root.classList.add('ui-panel');
    this.viewers.push(logViewer);
    this.root.append(ui.splitH([
      filters,
      logViewer.root,
    ]));
  }
}
