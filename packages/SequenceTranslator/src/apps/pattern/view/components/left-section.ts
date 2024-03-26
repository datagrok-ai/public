import * as ui from 'datagrok-api/ui';

import $ from 'cash-dom';
import {PatternDefaultsProvider} from '../../model/defaults-provider';
import {EventBus} from '../../model/event-bus';
import {PatternAppDataManager} from '../../model/external-data-manager';

import {PatternEditControlsManager} from './edit-block-controls';
import {PatternLoadControlsManager} from './load-block-controls';
import {TableControlsManager} from './bulk-convert/bulk-convert';

export class PatternAppLeftSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private defaults: PatternDefaultsProvider
  ) { };

  getLayout(): HTMLDivElement {
    const loadControlsManager = new PatternLoadControlsManager(
      this.eventBus,
      this.dataManager
    );

    const editControlsManager = new PatternEditControlsManager(
      this.eventBus,
      this.defaults
    );
    const tableControlsManager = new TableControlsManager(this.eventBus);

    const loadControls = loadControlsManager.createUIComponents();
    const editControls = editControlsManager.createUIComponents();
    const tableControls = tableControlsManager.createUIComponents();

    const loadControlsContainer = ui.div(loadControls);
    $(loadControlsContainer).css({'padding-bottom': '20px'});

    const form = ui.div([
      ...editControls,
      ...tableControls
    ], 'ui-form');

    const container = ui.div([
      loadControlsContainer,
      form
    ]);
    $(container).css({'padding': '25px'});

    const layout = ui.box(container, {style: {'maxWidth': '450px'}});
    return layout;
  }
}
