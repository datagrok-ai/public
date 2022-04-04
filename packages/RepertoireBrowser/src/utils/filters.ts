import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

export class CustomFilter {
  protected df: DG.DataFrame;
  protected ptmMap: StringMap;
  protected cdrMap: StringMap;
  protected referenceDf: DG.DataFrame;
  protected indexes: Int32Array;
  protected normalizedCDRMap: StringMap;

  constructor(
    df: DG.DataFrame,
    ptmMap: StringMap,
    cdrMap: StringMap,
    referenceDf: DG.DataFrame,
    indexes: Int32Array,
  ) {
    this.df = df;
    this.ptmMap = ptmMap;
    this.cdrMap = cdrMap;
    this.referenceDf = referenceDf,
    this.indexes = indexes;
  }

  //FIXME: can be preprocessed and saved
  protected normalizeCDRnames(v: string, ct: string): string {
    let cdrName = v.split('_')[0];
    cdrName = cdrName.charAt(0).toUpperCase() + cdrName.slice(1);
    this.normalizedCDRMap[`${ct} ${cdrName}`] = v;
    return cdrName;
  };

  create() {
    const ptmKeys = Object.keys(this.ptmMap);
    const ptmHKeys = ptmKeys.filter((v) => v.startsWith('H')).map((v) => v.slice(2));
    const ptmLKeys = ptmKeys.filter((v) => v.startsWith('L')).map((v) => v.slice(2));
    const cdrKeys = Object.keys(this.cdrMap);

    this.normalizedCDRMap = {};

    const cdrHKeys = cdrKeys.filter((v) => v.includes('CDRH')).map((v) => this.normalizeCDRnames(v, 'H'));
    const cdrLKeys = cdrKeys.filter((v) => v.includes('CDRL')).map((v) => this.normalizeCDRnames(v, 'L'));

    const options = this.options;
    const filterH = new FilterGroup(options, ptmHKeys, cdrHKeys, 'H').create();
    const filterL = new FilterGroup(options, ptmLKeys, cdrLKeys, 'L').create();

    filterH.style.margin = '10px';
    filterL.style.margin = '10px';
    return new DG.Widget(ui.divV([filterH, filterL]));
  }

  get options(): CustomFilterOptions {
    return {
      df: this.df,
      ptmMap: this.ptmMap,
      cdrMap: this.cdrMap,
      referenceDf: this.referenceDf,
      indexes: this.indexes,
      normalizedCDRMap: this.normalizedCDRMap,
    };
  }
}

class FilterGroup {
  protected probabilityHeader: HTMLDivElement;
  protected ptmInput: DG.InputBase;
  protected cdrInput: DG.InputBase;
  protected probabilityInput: DG.RangeSlider;

  protected options: CustomFilterOptions;

  protected ptmNames: string[];
  protected cdrNames: string[];
  protected chainType: string;

  constructor(
    options: CustomFilterOptions,
    ptmNames: string[],
    cdrNames: string[],
    chainType: string,
  ) {
    this.options = options;
    this.ptmNames = ptmNames;
    this.cdrNames = cdrNames;
    this.chainType = chainType;
  }

  protected onFilterApplyButton() {
    const cMin = parseFloat(this.probabilityInput.min.toFixed(3));
    const cMax = parseFloat(this.probabilityInput.max.toFixed(3));
    const ptmInputValue: string[] = this.ptmInput.value;
    const currentCdr = this.cdrInput.stringValue === 'None' ?
      'max' : this.options.cdrMap[this.options.normalizedCDRMap[`${this.chainType} ${this.cdrInput.stringValue}`]];

    this.options.df.filter.init((index: number) => {
      for (const chosenPTM of ptmInputValue) {
        const colName = this.options.ptmMap[`${this.chainType} ${chosenPTM}`];
        const binStr = this.options.referenceDf.get(colName, this.options.indexes[index]);

        if (typeof binStr === 'undefined' || binStr === '')
          return false;

        const cdrs = JSON.parse(binStr);
        if (!Object.keys(cdrs).includes(currentCdr))
          return false;

        const currentProbability = cdrs[currentCdr];
        if (typeof currentProbability === 'undefined' || currentProbability > cMax || currentProbability < cMin)
          return false;
      }

      return true;
    });
    this.options.df.filter.fireChanged();
  }

  protected onProbabilityInputValuesChanged() {
    const min = this.probabilityInput.min.toFixed(3);
    const max = this.probabilityInput.max.toFixed(3);
    this.probabilityHeader.textContent = `Probability: [${min}; ${max}]`;
  }

  create(): HTMLDivElement {
    this.probabilityHeader = ui.divText('Probability: [0.5; 1]');

    this.cdrNames.unshift('None');

    //@ts-ignore: method api is wrong
    this.ptmInput = ui.multiChoiceInput('PTM', [] as string[], this.ptmNames);
    this.cdrInput = ui.choiceInput('CDR', this.cdrNames[0], this.cdrNames);
    this.probabilityInput = ui.rangeSlider(0, 1, 0.5, 1);

    const probabilityHost = ui.divV([
      this.probabilityHeader,
      ui.divH([
        ui.inlineText(['0']),
        this.probabilityInput.root, ui.inlineText(['1']),
      ]),
    ]);

    this.probabilityInput.onValuesChanged.subscribe(this.onProbabilityInputValuesChanged.bind(this));

    $(this.probabilityInput.root).children('*').height('17px').css('max-width', '265px');

    const filterBtn = ui.button('Apply', this.onFilterApplyButton.bind(this));
    const heading = ui.h3(`${this.chainType} chain filter`);

    heading.style.margin = '0px';
    return ui.divV([heading, this.cdrInput.root, this.ptmInput.root, probabilityHost, filterBtn]);
  }
}

export type StringMap = {[key: string]: string};

type CustomFilterOptions = {
  df: DG.DataFrame;
  ptmMap: StringMap;
  cdrMap: StringMap;
  referenceDf: DG.DataFrame;
  indexes: Int32Array;
  normalizedCDRMap: StringMap;
};
