import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {_startInit} from './package';
import {MlbEvents} from './const';

type ChainTypeType = 'L' | 'H';

/**
 * Single-option categorical filter that demonstrates the concept of collaborative filtering:
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */
export class PtmFilter extends DG.Filter {
  cdrInput: DG.InputBase;

  chainType: ChainTypeType;
  ptmMap: { [key: string]: string; };
  cdrMap: { [key: string]: string; };
  refDf: DG.DataFrame;
  ptmKeys: string[];
  cdrKeys: string[];
  normalizedCDRMap: { [key: string]: string };
  indexes: Int32Array;
  cMin: number;
  cMax: number;
  ptmInputValue: string[];
  currentCdr: string;

  constructor(ptmMap: { [key: string]: string }, cdrMap: { [key: string]: string }, refDf: DG.DataFrame) {
    super();
    this.root = ui.divV(null, 'd4-mlb-filter');
    this.subs = [];
    this.ptmMap = ptmMap;
    this.cdrMap = cdrMap;
    this.refDf = refDf;

    this.ptmKeys = [...new Set(Object.keys(ptmMap).map((v) => v.slice(2)))];
    this.normalizedCDRMap = {};
    //FIXME: can be preprocessed and saved
    const normalizeCDRnames = (v: string) => {
      let cdrName = v.split('_')[0];
      cdrName = cdrName.charAt(0).toUpperCase() + cdrName.slice(1);
      this.normalizedCDRMap[`L ${cdrName}`] = v;
      this.normalizedCDRMap[`H ${cdrName}`] = v;
      return cdrName;
    };

    this.cdrKeys = [...new Set(Object.keys(this.cdrMap).map((v) => normalizeCDRnames(v)))];
    this.cdrKeys.unshift('None');

    this.cMin = 0;
    this.cMax = 1;
    this.ptmInputValue = [];
    this.currentCdr = this.cdrKeys[0];
    this.chainType = 'L';
  }

  get caption() {
    return 'PTM filter';
  }

  get isFiltering() {
    return true;
  }

  get filterSummary() {
    return `${this.ptmInputValue} in ${this.currentCdr} with probability in [${this.cMin}, ${this.cMax}]`;
  }

  attach(dataFrame: DG.DataFrame) {
    try {
      console.debug(`MLB: PtmFilter.attach() start, ${((Date.now() - _startInit) / 1000).toString()}`);
      super.attach(dataFrame);

      const tempDf = this.refDf.clone(null, ['v_id']);
      tempDf.columns.addNewInt('index').init((i) => i);
      this.indexes = (dataFrame.clone(null, ['v id']).join(tempDf, ['v id'], ['v_id'], [], ['index'], 'left', false)
        .getCol('index').getRawData() as Int32Array);

      this.render();
      console.debug(`MLB: PtmFilter.attach() end, ${((Date.now() - _startInit) / 1000).toString()}`);
    } catch (err: unknown) {
      console.error(err instanceof Error ? err.message : (err as Object).toString());
    }
  }

  applyState(state: any) {
    super.applyState(state);

    // OnPropertyChanged() is not called when setting a property in TableView.filters() argument,
    // but applyState() is called with that value
    if ('currentCdr' in state) {
      let value: string = state['currentCdr'];
      const key: string = this.cdrKeys.find((v) => v.toUpperCase() == value.toUpperCase());
      console.debug(`MLB: PtmFilter.applyState(${MlbEvents.CdrChanged}) value ="${value}" -> key="${key}".`);
      this.currentCdr = key ?? 'None';
    }

    this.render();
  }

  detach() {
    super.detach();
    console.log('MLB: PtmFilter.detach() filter detached');
  }

  applyFilter() {
    try {
      console.debug(`MLB: PtmFilter.applyFilter() start, ${((Date.now() - _startInit) / 1000).toString()}`);

      const getStateFor = (index: number) => {
        for (const chosenPTM of this.ptmInputValue) {
          const binStr = this.refDf.get(this.ptmMap[`${this.chainType} ${chosenPTM}`], this.indexes[index]);
          if (typeof binStr === 'undefined' || binStr === '')
            return false;

          const cdrs = JSON.parse(binStr);
          if (!Object.keys(cdrs).includes(this.currentCdr))
            return false;

          const currentProbability = cdrs[this.currentCdr];
          if (
            typeof currentProbability === 'undefined' ||
            currentProbability > this.cMax || currentProbability < this.cMin
          ) return false;
        }

        return true;
      };

      const filter = this.dataFrame.filter;
      const rowCount = this.dataFrame.rowCount;

      for (let i = 0; i < rowCount; i++)
        filter.set(i, filter.get(i) && getStateFor(i), false);

      filter.fireChanged();
      console.debug(`MLB: PtmFilter.applyFilter() end, ${((Date.now() - _startInit) / 1000).toString()}`);
    } catch (err: unknown) {
      console.error(err instanceof Error ? err.message : (err as Object).toString());
    }
  }

  render() {
    console.debug(`MLB: PtmFilter.render() start, ${((Date.now() - _startInit) / 1000).toString()}`);

    $(this.root).empty();

    const chainTypeInput = ui.choiceInput('Chain type', 'L', ['L', 'H'], () => {
      this.chainType = chainTypeInput.stringValue as ChainTypeType;
      this.dataFrame.rows.requestFilter();
    });
    const ptmInputLabel = ui.label('PTM');
    ptmInputLabel.style.marginTop = '5px';
    //@ts-ignore: method api is wrong
    const ptmInput = ui.multiChoiceInput('PTM', this.ptmInputValue, this.ptmKeys, () => {
      this.ptmInputValue = ptmInput.value;
      this.dataFrame.rows.requestFilter();
    });
    ptmInput.root.style.width = 'max-content';
    const ptmInputCash = $(ptmInput.root);
    ptmInputCash.find('.ui-input-multi-choice-checks > div > input.ui-input-editor[type="checkbox"]')
      .each((_i, element) => {
        element.style.margin = '0';
      });
    ptmInputCash.find('.ui-input-multi-choice-checks > div').each((_i, element) => {
      ui.tooltip.bind(element, `${element.lastChild.textContent}`);
    });
    ptmInputCash.children().first().remove();
    this.cdrInput = ui.choiceInput('CDR', this.currentCdr, this.cdrKeys, () => {
      this.currentCdr = this.cdrInput.stringValue === 'None' ?
        'max' : this.cdrMap[this.normalizedCDRMap[`${this.chainType} ${this.cdrInput.stringValue}`]];
      this.dataFrame.rows.requestFilter();
    });

    const probabilityHeader = ui.divText(this.getProbabilityText());
    probabilityHeader.style.marginTop = '5px';
    probabilityHeader.style.marginBottom = '5px';
    const probabilityInput = ui.rangeSlider(0, 1, this.cMin, this.cMax, false, 'thin_barbell');

    const rsMin = ui.inlineText(['0']);
    rsMin.style.marginRight = '5px';
    const rsMax = ui.inlineText(['1']);
    rsMax.style.marginLeft = '5px';
    const probabilityHost = ui.divV([
      probabilityHeader, ui.divH([rsMin, probabilityInput.root, rsMax])]);

    probabilityInput.onValuesChanged.subscribe((_) => {
      this.cMin = parseFloat(probabilityInput.min.toFixed(3));
      this.cMax = parseFloat(probabilityInput.max.toFixed(3));
      probabilityHeader.textContent = this.getProbabilityText();
      this.dataFrame.rows.requestFilter();
    });
    $(probabilityInput.root).children('*').height('17px').css('max-width', '265px');

    this.root = ui.divV([chainTypeInput.root, this.cdrInput.root, ptmInputLabel, ptmInput.root, probabilityHost]);
    this.root.style.margin = '10px';

    console.debug(`MLB: PtmFilter.render() end, ${((Date.now() - _startInit) / 1000).toString()}`);
  }

  private getProbabilityText(): string {
    return `Probability: [${this.cMin}; ${this.cMax}]`;
  }
}
