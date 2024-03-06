/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { PatternConfiguration, StrandType } from '../../model/types';
import {PatternConfigurationManager} from '../../model/pattern-config-manager';
import {EventBus} from '../../model/event-bus';
import {NucleotidePatternSVGRenderer} from './svg-renderer';

export class SvgDisplayManager {
  private svgDisplayDiv = ui.div([]);

  private constructor(
    private eventBus: EventBus
  ) {
    eventBus.patternStateChanged$.subscribe(() => this.updateSvgContainer());
  }

  static createSvgDiv(eventBus: EventBus): HTMLDivElement {
    const displayManager = new SvgDisplayManager(eventBus);
    return displayManager.svgDisplayDiv;
  }

  private updateSvgContainer(): void {
    $(this.svgDisplayDiv).empty();
    const patternConfig = PatternConfigurationManager.getConfig(this.eventBus);
    const image = this.createSvg(patternConfig);
    this.svgDisplayDiv.append(image);
  }

  private createSvg(patternConfig: PatternConfiguration) {
    const renderer = new NucleotidePatternSVGRenderer(patternConfig);
    const svg = renderer.renderPattern();
    return svg;
  }
}
