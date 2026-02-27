import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {UaView} from './ua';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {UaToolbox} from '../ua-toolbox';

const filters = ui.box();
filters.style.maxWidth = '250px';

const users: Array<any> = [];
// const tags = ui.divH([]);

export class LogView extends UaView {
  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Log';
  }

  async initViewers(path?: string): Promise<void> {
    grok.dapi.users.list().then(
      (user) => user.map((u) => {
        const d = {
          avatar: u.picture,
          name: u.friendlyName,
          data: u,
        };
        users.push(d);
      //add object of users with name and picture element;
      }));

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
        viewer.col('user')!.cellType = 'html';
        viewer.col('user')!.width = 30;
        viewer.col('source')!.width = 30;
        viewer.col('ugid')!.visible = false;

        viewer.onCellPrepare(async function(gc) {
          if (gc.gridColumn.name === 'event_description')
            gc.style.font = '13px monospace';

          if (gc.gridColumn.name === 'event_time' || gc.gridColumn.name === 'id') {
            gc.style.textColor = 0xFFB8BAC0;
            gc.style.font = '13px Roboto';
          }

          if (gc.gridColumn.name === 'user') {
            for (let i = 0; i < users.length; i++) {
              if (users[i].name == gc.cell.value) {
                const img = ui.div();
                img.style.width = '20px';
                img.style.height = '20px';
                img.style.backgroundSize = 'contain';
                img.style.margin = '5px 0 0 5px';
                img.style.borderRadius = '100%';
                if (gc.cell.value != 'Test')
                  img.style.backgroundImage = users[i].avatar.style.backgroundImage;
                else
                  img.style.backgroundImage = 'url(/images/entities/grok.png);';
                img.addEventListener('click', ()=>{
                  grok.shell.o = users[i].data;
                });
                gc.style.element = ui.tooltip.bind(img, users[i].name);
              }
            }
          }
        });

        viewer.onCellRender.subscribe(function(args) {
          if (args.cell.gridColumn.name == 'source') {
            args.g.beginPath();
            switch (args.cell.cell.value) {
            case 'error': args.g.fillStyle = '#d62727';
              break;
            case 'audit': args.g.fillStyle = '#1f77b4';
              break;
            case 'data-query': args.g.fillStyle = '#9467bd';
              break;
            case 'function': args.g.fillStyle = '#2ba02b';
              break;
            case 'function-package': args.g.fillStyle = '#97df8a';
              break;
            case 'script': args.g.fillStyle = '#ffbb78';
              break;
            default: args.g.fillStyle = '#7f7f7f';
              break;
            }
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

const filtersStyle = {
  columnNames: ['source', 'user', 'event_time'],
};
