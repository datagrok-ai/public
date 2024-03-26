import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {SvgDisplayManager} from '../svg-utils/svg-display-manager';

import {EventBus} from '../../model/event-bus';
import {PatternAppDataManager} from '../../model/external-data-manager';
import {isOverhangNucleotide} from '../../model/utils';
import {PatternNameExistsError, PatternExistsError} from '../../model/types';

import $ from 'cash-dom';
import {BooleanInput} from '../types';

export class PatternAppRightSection {
  private svgDisplay: HTMLDivElement;

  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager
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
    this.dataManager.savePatternToUserStorage()
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

class NumericLabelVisibilityControls {
  private togglesContainer: HTMLDivElement = ui.div([]);

  constructor(
    private eventBus: EventBus
  ) {
    this.eventBus.uniqueNucleotidesChanged$().subscribe(() => {
      this.updateContainer();
    });
  }

  getContainer(): HTMLDivElement {
    return this.togglesContainer;
  }

  private updateContainer(): void {
    const newInputs = this.createInputs();
    $(this.togglesContainer).empty();
    $(this.togglesContainer).append(newInputs);
  }

  private createInputs(): HTMLDivElement {
    const uniqueNucleotideBases = this.eventBus.getUniqueNucleotides();
    const nucleotidesWithoutOverhangs = uniqueNucleotideBases.filter((n) => !isOverhangNucleotide(n));

    const inputBases = nucleotidesWithoutOverhangs.map(
      (nucleotide: string) => this.createSingleInput(nucleotide)
    );

    inputBases.sort(
      (inputA, inputB) => inputA.captionLabel.textContent!.localeCompare(inputB.captionLabel.textContent!)
    );

    return ui.divH(inputBases.map((input) => input.root));
  }

  private createSingleInput(nucleotide: string): BooleanInput {
    const initialValue = this.eventBus.getModificationsWithNumericLabels().includes(nucleotide);
    const input = ui.boolInput(
      nucleotide,
      initialValue,
      (value: boolean) => this.handleNumericLabelToggle(nucleotide, value)
    );
    $(input.root).css('padding-right', '20px');

    input.setTooltip(`Show numeric labels for ${nucleotide}`);

    return input;
  }

  private handleNumericLabelToggle(nucleotide: string, isVisible: boolean): void {
    const labelledNucleotides = this.eventBus.getModificationsWithNumericLabels();
    const hasNumericLabel = labelledNucleotides.includes(nucleotide);
    if (hasNumericLabel === isVisible)
      return;

    const newArray = isVisible ? labelledNucleotides.concat(nucleotide) :
      labelledNucleotides.filter((n) => n !== nucleotide);
    this.eventBus.updateModificationsWithNumericLabels(newArray);
  }
}
