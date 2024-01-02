/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, StrandType, OTHER_USERS } from '../model/const';

import {StringInput, NumberInput} from './types';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternConfigurationManager} from '../model/pattern-state-manager';
import * as rxjs from 'rxjs';
// WARNING: for some reason, cannot use rxjs.operators.debounceTime, although
// webpack.config.js is configured to use rxjs.operators as rxjs.operators
import $ from 'cash-dom';

export class PatternAppRightSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    private patternConfiguration: PatternConfigurationManager,
  ) { };

  getLayout(): HTMLDivElement {
    const svgDisplay = new SvgDisplayManager().createUI();
    const layout = ui.panel([
      svgDisplay,
      // numericLabelTogglesContainer,
      // generateDownloadAndEditControls(),
      // generateStrandSectionDisplays(),
      // ui.h1('Additional modifications'),
      // ui.form([
      //   terminalModificationInputs[SENSE_STRAND][FIVE_PRIME_END],
      //   terminalModificationInputs[SENSE_STRAND][THREE_PRIME_END],
      // ]),
      // asModificationDiv,
      ], {style: {overflowX: 'scroll', padding: '12px 24px'}});
    return layout;
  }

  private createNumericLabelTogglesContainer(): HTMLElement {
    const defaultNucleotideBase = this.patternConfiguration.getDefaultNucleotideBase();
    const numericLabelTogglesContainer = ui.divH([
      ui.boolInput(defaultNucleotideBase, true, (value: boolean) => {
        updateListOfModificationsWithNumericLabels(value);
        refreshSvgDisplay();
        refreshOutputExamples();
      }).root,
    ]);
    return numericLabelTogglesContainer;
  }
}

class SvgDisplayManager {
  constructor() { }

  private svgDisplayDiv = ui.div([]);
  createUI(): HTMLElement {
    return this.svgDisplayDiv;
  }
}
