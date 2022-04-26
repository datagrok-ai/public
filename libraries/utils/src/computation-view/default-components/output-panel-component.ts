/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Component, ComputationViewStateService} from '../common/service-interfaces';
import {DefaultComputationViewStateService} from '../default-services/state-service';

export class DefaultOutputPanelComponent implements Component {
  public root: HTMLElement = ui.div([]);

  constructor(
    private state: ComputationViewStateService = new DefaultComputationViewStateService,
  ) {
    this.state.onComputationCompleted.subscribe(() => this.render());
  }

  render(): void {
    const outputs = this.state.lastCall?.outputs as Record<string, any>;
    console.log(Object.entries(outputs));
    const panel = ui.accordion('Output data');
    panel.addPane('Output data', () => {
      return ui.divV(
        Object.entries(outputs).map(([key, val]) => ui.span([`${key}: `, `${val}`])),
      );
    });
    this.root.replaceWith(panel.root);
  }
}
