/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import _ from 'lodash';
import $ from 'cash-dom';

import {EventBus} from '../model/event-bus';
import {AXOLABS_STYLE_MAP} from '../../common/data-loader/json-loader';
import {BooleanInput, StringInput} from './types';
import {PatternConfiguration, PhosphorothioateLinkageFlags, StrandType} from '../model/types';
import {STRAND, STRANDS, STRAND_LABEL} from '../model/const';

export let isDialogOpen = false;

export class PatternEditorDialog {
  private initialPatternConfig: PatternConfiguration;

  constructor(
    private eventBus: EventBus,
  ) {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
  }

  open(): void {
    isDialogOpen = true;
    this.createDialog().show();
  }

  private createDialog(): DG.Dialog {
    const editorBody = ui.divV([]);
    this.eventBus.updatePatternEditor$.subscribe(() => {
      this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
      const header = new HeaderControls(this.eventBus, this.initialPatternConfig).getPhosphorothioateLinkageControls();
      const controls = new StrandControls(this.eventBus, this.initialPatternConfig).create();

      $(editorBody).empty();
      $(editorBody).append(header, controls);
    });

    const dialog = ui.dialog('Edit pattern')
    .add(editorBody)
    .onOK(() => {})
    .onCancel(() => this.resetToInitialState());

    dialog.onClose.subscribe(() => isDialogOpen = false);

    return dialog;
  }

  private resetToInitialState(): void {
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

class HeaderControls {
  constructor(
    private eventBus: EventBus,
    private initialPatternConfig: PatternConfiguration
  ) { }

  getPhosphorothioateLinkageControls(): HTMLDivElement {
    const container = ui.divV([
      ui.h1('PTO'),
      ui.divH([
        this.createAllPtoActivationInput().root,
        ...this.createFirstPtoInputs().map((input) => input.root),
      ], {style:{gap:'12px'}})
    ])

    return container;
  }

  private areAllPtoLinkagesSet(flags: PhosphorothioateLinkageFlags): boolean {
    const totalNumberOfPTOFlags = STRANDS.map(
      (strand) => flags[strand].filter((flag) => flag).length
    )
      .reduce((a, b) => a + b, 0);
    const totalNumberOfNucleotides = STRANDS.map((strand) => this.initialPatternConfig.nucleotideSequences[strand].length).reduce((a, b) => a + b, 0);

    // WARNING: +2 because there are +1 more PTO flags in each strand than there are nucleotides
    return totalNumberOfPTOFlags === totalNumberOfNucleotides + 2;
  }

  private createAllPtoActivationInput(): BooleanInput {
    const flags = this.initialPatternConfig.phosphorothioateLinkageFlags;
    const initialValue = this.areAllPtoLinkagesSet(flags);
    const allPtoActivationInput = ui.boolInput('All PTO', initialValue);

    allPtoActivationInput.onInput(() => {
      const value = allPtoActivationInput.value!;
      this.eventBus.setAllPTOLinkages(value);
    });

    this.eventBus.phosphorothioateLingeFlagsChanged$.subscribe(() => {
      const flags = this.eventBus.getPatternConfig().phosphorothioateLinkageFlags;
      const newValue = this.areAllPtoLinkagesSet(flags);
      allPtoActivationInput.value = newValue;
    });

    this.addStyleToPtoInput(allPtoActivationInput);
    ui.tooltip.bind(allPtoActivationInput.captionLabel, 'Include all phosphothioates in the pattern');

    return allPtoActivationInput;
  }

  private addStyleToPtoInput(allPtoActivationInput: BooleanInput): void {
    const label = allPtoActivationInput.captionLabel;
    label.classList.add('ui-label-right');
    Object.assign(label.style, {
      textAlign: 'left',
      maxWidth: '100px',
      minWidth: '40px',
      width: 'auto'
    });
  }

  private createFirstPtoInputs(): BooleanInput[] {
    return STRANDS.map((strand) => {
      const initialValue = this.isFirstPtoActive(strand);
      const firstPtoInput = ui.boolInput(`First ${strand} PTO`, initialValue);

      firstPtoInput.onInput(() => {
        const value = firstPtoInput.value!;
        this.eventBus.setPhosphorothioateLinkageFlag(strand, 0, value);
      });

      this.eventBus.phosphorothioateLingeFlagsChanged$.subscribe((flags) => {
        const newValue = flags[strand][0];
        firstPtoInput.value = newValue;
      });
      
      this.addStyleToPtoInput(firstPtoInput);
      ui.tooltip.bind(firstPtoInput.captionLabel, `Include first phosphothioate in ${strand}`);
      return firstPtoInput;
    });
  }

  private isFirstPtoActive(strand: STRAND): boolean {
    return this.initialPatternConfig.phosphorothioateLinkageFlags[strand][0];
  }
}

class StrandControls {
  private displayedInputLabels: Map<StrandType, string[]>;

  constructor(
    private eventBus: EventBus,
    private initialPatternConfig: PatternConfiguration,
  ) {
    this.eventBus.nucleotideSequencesChanged$.subscribe(() => {
      this.displayedInputLabels = this.computeDisplayedInputLabels();
    });
  }

  create(): HTMLDivElement {
    const inputPanels = STRANDS.map((strand) => this.constructControlsPanel(strand));

    const container = ui.divH(inputPanels, {style:{gap:'24px'}})
    return container;
  }

  private constructControlsPanel(strand: StrandType): HTMLDivElement {
    const header = this.constructHeader();
    const modificationControls = this.createControls(strand);

    const container = ui.block([
      ui.h1(`${STRAND_LABEL[strand]}`),
      header,
      modificationControls,
    ], {style: {paddingTop: '12px'}});

    return container;
  }

  private constructHeader() {
    return ui.divH([
      ui.div([ui.divText('#')], {style: {width: '20px'}}),
      ui.block75([ui.divText('Modification')]),
      ui.div([ui.divText('PTO')]),
    ]);
  }

  private createControls(strand: StrandType): HTMLDivElement {
    const nucleobaseInputs = this.createNucleobaseInputs(strand);
    const labels = this.createLabelDivs(strand);
    const ptoLinkageInputs = this.createPTOFlagInputs(strand);
    const container = ui.div(nucleobaseInputs.map(
      (nucleobaseInput, idx) => {
        return ui.divH([
          labels[idx],
          nucleobaseInput.root,
          ptoLinkageInputs[idx].root,
        ], {style: {alignItems: 'center'}});
    }));
    return container;
  }

  private createNucleobaseInputs(strand: StrandType): StringInput[] {
    const nucleotideBaseChoices: string[] = Object.keys(AXOLABS_STYLE_MAP).sort();
    const nucleotides = this.eventBus.getNucleotideSequences()[strand];
    const choiceInputs = nucleotides.map((nucleotide, index) => {
      const input = ui.choiceInput<string>('', nucleotide, nucleotideBaseChoices);
      input.onInput(() => {
        const newValue = input.value!;
        this.eventBus.setNucleotideBase(strand, index, newValue);
      });
      return input;
    });

    return choiceInputs;
  }
  
  private createPTOFlagInputs(strand: StrandType): BooleanInput[] {
    const ptoLinkageFlags = this.eventBus.getPhosphorothioateLinkageFlags()[strand].slice(1);
    const ptoLinkageInputs = ptoLinkageFlags.map((flag, index) => {
      const input = ui.boolInput('', flag);
      input.onInput(() => {
        const newValue = input.value!;
        this.eventBus.setPhosphorothioateLinkageFlag(strand, index + 1, newValue);
      })
      return input;
    });

    return ptoLinkageInputs;
  }

  private computeDisplayedInputLabels(): Map<StrandType, string[]> {
    const nucleotides = this.eventBus.getNucleotideSequences();
    const labels = new Map<StrandType, string[]>();
    STRANDS.forEach((strand) => {
      let counter = 1;
      const strandNucleotides = nucleotides[strand];
      const strandLabels = strandNucleotides.map((nucleotide) => {
        if (this.isOverhangNucleotide(nucleotide)) {
          return '';
        }
        const label = String(counter);
        counter++;
        return label;
      });
      labels.set(strand, strandLabels);
    });

    return labels;
  }

  private createLabelDivs(strand: StrandType): HTMLElement[] {
    const labels = this.createLabels(strand);
    const labelDivs = labels.map((label) => ui.div([label], {style: {width: '20px'}}));

    this.eventBus.nucleotideSequencesChanged$.subscribe(() => {
      const newLabels = this.createLabels(strand);
      newLabels.forEach((newLabel, index) => {
        $(labelDivs[index]).empty();
        $(labelDivs[index]).append(newLabel);
      });
    });

    return labelDivs;
  }

  private createLabels(strand: StrandType): HTMLLabelElement[] {
    const nucleotides = this.eventBus.getNucleotideSequences()[strand];
    const labels = nucleotides.map((_, index) => {
      const labelText = this.displayedInputLabels.get(strand)![index];
      return ui.label(labelText);
    });

    return labels;
  }


  private isOverhangNucleotide(modification: string): boolean {
    return modification.endsWith('(o)');
  }

  // private generateSingleModificationControlGroup(strand: string, nucleotideCounter: number, index: number) {
  //   const labelText = isOverhangNucleotide(getBaseInputValue(strand, index)) ? '' : String(nucleotideCounter);
  //   const labelUI = ui.div([ui.label(labelText)], {style: {width: '20px'}});
  //   const baseInputUI = ui.block75([nucleobaseInputs[strand][index].root]);
  //   const ptoLinkageUI = ui.div([ptoLinkageInputs[strand][index]]);

  //   return ui.divH([labelUI, baseInputUI, ptoLinkageUI], {style: {alignItems: 'center'}});
  // }
}
