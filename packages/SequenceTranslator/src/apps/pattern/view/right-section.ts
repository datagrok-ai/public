/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { STRAND, STRAND_LABEL, STRANDS, OTHER_USERS } from '../model/const';
import { StrandType } from '../model/types';

import {StringInput, NumberInput} from './types';
import {NucleotidePatternSVGRenderer} from '../view/svg-utils/svg-renderer';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import * as rxjs from 'rxjs';
import $ from 'cash-dom';

export class PatternAppRightSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
  ) { };

  getLayout(): HTMLDivElement {
    const svgDisplay = new SvgDisplayManager().create();
    const layout = ui.panel([
      svgDisplay,
      // numericLabelTogglesContainer,
      // generateDownloadAndEditControls(),
      // generateStrandSectionDisplays(),
      // ui.h1('Additional modifications'),
      // ui.form([
      //   terminalModificationInputs[STRAND.SENSE][FIVE_PRIME_END],
      //   terminalModificationInputs[STRAND.SENSE][THREE_PRIME_END],
      // ]),
      // asModificationDiv,
      ], {style: {overflowX: 'scroll', padding: '12px 24px'}});
    return layout;
  }

  // private createNumericLabelTogglesContainer(): HTMLElement {
  //   const defaultNucleotideBase = this.patternConfiguration.getDefaultNucleotideBase();
  //   const numericLabelTogglesContainer = ui.divH([
  //     ui.boolInput(defaultNucleotideBase, true, (value: boolean) => {
  //       updateListOfModificationsWithNumericLabels(value);
  //       refreshSvgDisplay();
  //       refreshOutputExamples();
  //     }).root,
  //   ]);
  //   return numericLabelTogglesContainer;
  // }
}

class SvgDisplayManager {
  private svgDisplayDiv = ui.div([]);
  // private renderer = new NucleotidePatternSVGRenderer(patternConfiguration);
  // const svg = renderer.renderPattern();

  constructor() { }

  create(): HTMLElement {
    return this.svgDisplayDiv;
  }
}
