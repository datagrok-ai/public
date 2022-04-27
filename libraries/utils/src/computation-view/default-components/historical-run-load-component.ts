/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Component, HistoricalRunService} from '../common/service-interfaces';
import {tokens} from '../common/inject-tokens';

export class DefaultHistoricalRunComponent implements Component {
  public root: HTMLElement = ui.div();

  public static inject = [
    tokens.historyService,
  ] as const;

  constructor(private historicalRunService: HistoricalRunService) {
    this.render();
  }

  render(): void {
    this.root = ui.iconFA('history', () => {
      this.historicalRunService.pullRuns().then(async (historicalRuns) => {
        const menu = DG.Menu.popup();
        let i = 0;
        for (const run of historicalRuns) {
          //@ts-ignore
          menu.item(`${run.func.friendlyName} — ${i}`, async () => this.historicalRunService.loadRun(run.id));
          i++;
        }
        menu.show();
      });
    });
  }
}
