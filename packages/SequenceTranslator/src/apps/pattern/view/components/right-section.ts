import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {EventBus} from '../../model/event-bus';
import {DataManager} from '../../model/data-manager';
import {PatternExistsError, PatternNameExistsError} from '../../model/types';
import {SvgDisplayManager} from '../svg-utils/svg-display-manager';
import {NumericLabelVisibilityControls} from './numeric-label-visibility-controls';

export class PatternAppRightSection {
  private svgDisplay: HTMLDivElement;

  constructor(
    private eventBus: EventBus,
    private dataManager: DataManager
  ) {
    this.svgDisplay = SvgDisplayManager.createSvgDiv(eventBus);
  };

  getLayout(): HTMLDivElement {
    const numericLabelTogglesContainer = new NumericLabelVisibilityControls(this.eventBus).getContainer();
    const downloadAndEditControls = this.generateDownloadControls();
    const layout = ui.panel([
      this.svgDisplay,
      numericLabelTogglesContainer,
      downloadAndEditControls,
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

  private generateDownloadControls(): HTMLDivElement {
    return ui.divH([
      this.createSavePatternButton(),
      this.createDownloadPngButton(),
    ], {style: {gap: '12px', marginTop: '12px'}});
  }

  private createDownloadPngButton(): HTMLButtonElement {
    const svgDownloadButton = ui.button('Get PNG', () => this.eventBus.requestSvgSave());

    ui.tooltip.bind(svgDownloadButton, 'Download pattern as PNG');

    return svgDownloadButton;
  }

  private createSavePatternButton(): HTMLButtonElement {
    const savePatternButton = ui.button('Save pattern', () => this.processSaveButtonClick());

    ui.tooltip.bind(savePatternButton, 'Save pattern to user storage');

    return savePatternButton;
  }

  private processSaveButtonClick(): void {
    const patternName = this.eventBus.getPatternName();
    if (patternName === '') {
      grok.shell.warning(`Insert pattern name`);
      return;
    }
    this.dataManager.savePatternToUserStorage(this.eventBus)
      .then(() => {
        grok.shell.info(`Pattern ${patternName} saved`);
      })
      .catch((e) => {
        if (e instanceof PatternNameExistsError) {
          new OverwritePatternDialog(this.eventBus).show();
        } else if (e instanceof PatternExistsError) {
          grok.shell.warning(ui.div([
            ui.divText(`Pattern already exists`),
            ui.button('Load', () => {}),
          ]));
        } else {
          console.error('Error while saving pattern', e);
        }
      });
  }
}

class OverwritePatternDialog {
  constructor(
    private eventBus: EventBus
  ) { }

  show(): void {
    const patternName = this.eventBus.getPatternName();
    const dialog = ui.dialog(`Pattern ${patternName} already exists`);
    dialog.add(ui.divText(
      `Pattern ${patternName} already exists. Do you want to overwrite it?`
    ));
    dialog.show();

    dialog.onOK(() => this.processOverwriteNamesakePattern());
  }

  private processOverwriteNamesakePattern(): void {
    const patternName = this.eventBus.getPatternName();
    grok.shell.info(`Pattern ${patternName} overwritten`);
    // this.dataManager.overwritePatternInUserStorage()
    //   .then(() => {
    //     grok.shell.info(`Pattern ${patternName} overwritten`);
    //   })
    //   .catch((e) => {
    //     console.error('Error while overwriting pattern', e);
    //   });
  }
}

