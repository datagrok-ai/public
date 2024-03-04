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
  private configManager: PatternConfigurationManager;

  private constructor(
    eventBus: EventBus,
  ) {
    this.configManager = new PatternConfigurationManager(eventBus);
    eventBus.patternStateChanged$.subscribe(() => this.updateSvgContainer());
  }

  static createSvgDiv(eventBus: EventBus): HTMLElement {
    const manager = new SvgDisplayManager(eventBus);
    manager.updateSvgContainer();
    return manager.svgDisplayDiv;
  }

  private updateSvgContainer(): void {
    $(this.svgDisplayDiv).empty();
    const patternConfig = this.configManager.getConfig();
    const image = this.createSvg(patternConfig);
    this.svgDisplayDiv.append(image);
  }

  private createSvg(patternConfig: PatternConfiguration) {
    const renderer = new NucleotidePatternSVGRenderer(patternConfig);
    const svg = renderer.renderPattern();
    return svg;
  }
}
