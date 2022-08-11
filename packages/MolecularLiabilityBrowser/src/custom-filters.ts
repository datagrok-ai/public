import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';

import {_startInit} from './package';
import {MlbEvents} from './const';
import {Subscription, Unsubscribable} from 'rxjs';

type ChainTypeType = 'L' | 'H';

type PtmMapType = { [key: string]: string };
type CdrMapType = { [key: string]: string };

/**
 * Single-option categorical filter that demonstrates the concept of collaborative filtering:
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */
export class PtmFilter extends DG.Filter {
  predictedCdrInput: DG.InputBase;

  chainType: ChainTypeType;
  predictedPtmMap: PtmMapType;
  predictedCdrMap: CdrMapType;
  observedPtmMap: PtmMapType;
  observedCdrMap: CdrMapType;

  predictedPtmDf: DG.DataFrame;
  observedPtmDf: DG.DataFrame;

  predictedPtmKeys: string[];
  observedPtmKeys: string[];

  predictedCdrKeys: string[];
  observedCdrKeys: string[];

  normalizedCDRMap: { [key: string]: string };

  predictedIndexes: Int32Array;
  observedIndexes: Int32Array;

  predictedProbabilityMin: number;
  predictedProbabilityMax: number;
  observedProbabilityMin: number;
  observedProbabilityMax: number;

  // TODO: PTM checks to use symbols
  predictedPtmInputValue: string[];
  observedPtmInputValue: string[];

  currentCdr: string;

  viewSubs: Unsubscribable[];

  constructor(
    predictedPtmMap: PtmMapType, predictedCdrMap: CdrMapType,
    observedPtmMap: PtmMapType, observedCdrMap: CdrMapType
  ) {
    super();
    this.root = ui.divV(null, 'd4-mlb-filter');
    this.subs = [];
    this.predictedPtmMap = predictedPtmMap;
    this.predictedCdrMap = predictedCdrMap;
    this.observedPtmMap = observedPtmMap;
    this.observedCdrMap = observedCdrMap;

    this.predictedPtmKeys = [...new Set(Object.keys(predictedPtmMap).map((v) => v.slice(2)))];
    this.observedPtmKeys = [...new Set(Object.keys(observedPtmMap).map((v) => v.slice(2)))];
    this.normalizedCDRMap = {};
    //FIXME: can be preprocessed and saved
    const normalizeCDRnames = (v: string) => {
      let cdrName = v.split('_')[0];
      cdrName = cdrName.charAt(0).toUpperCase() + cdrName.slice(1);
      this.normalizedCDRMap[`L ${cdrName}`] = v;
      this.normalizedCDRMap[`H ${cdrName}`] = v;
      return cdrName;
    };

    this.predictedCdrKeys = [...new Set(Object.keys(this.predictedCdrMap).map((v) => normalizeCDRnames(v)))];
    this.predictedCdrKeys.unshift('None');

    this.observedCdrKeys = [...new Set(Object.keys(this.observedCdrMap).map((v) => normalizeCDRnames(v)))];
    this.observedCdrKeys.unshift('None');

    this.predictedProbabilityMin = 0;
    this.predictedProbabilityMax = 1;
    this.observedProbabilityMin = 0;
    this.observedProbabilityMax = 1;

    this.predictedPtmInputValue = [];
    this.observedPtmInputValue = [];

    this.currentCdr = this.predictedCdrKeys[0];
    this.chainType = 'L';
  }

  get caption() {
    return 'PTM filter';
  }

  get isFiltering() {
    return true;
  }

  get filterSummary() {
    return `${this.predictedPtmInputValue} in ${this.currentCdr} with probability ` +
      `in [${this.predictedProbabilityMin}, ${this.predictedProbabilityMax}]`;
  }

  /** some checks to investigate data structure and ensure key consistency
   * @param {PtmMapType} ptmMap
   * @param {CdrMapType} cdrMap
   * @param {DG.DataFrame} refDf Data frame to check for consistency with maps
   */
  static checkRefDf(ptmMap: PtmMapType, cdrMap: CdrMapType, refDf: DG.DataFrame) {
    const cdrSymbolSet = new Set(Object.values(cdrMap));
    for (const [ptmName, ptmSymbol] of Object.entries(ptmMap)) {
      for (let rowI = 0; rowI < refDf.rowCount; rowI++) {
        const refTxt = refDf.get(ptmSymbol, rowI);
        if (refTxt) {
          const refVal: { [cdrSymbol: string]: number } = JSON.parse(refTxt);

          if (!('max' in refVal))
            throw new Error('There is no "max".');

          // const max: number = refVal['max'];
          delete refVal['max'];

          for (const cdrSymbol of Object.keys(refVal)) {
            if (!cdrSymbolSet.has(cdrSymbol))
              throw new Error(`Unexpected cdrSymbol="${cdrSymbol}".`);
          }
        }
      }
    }
  }

  attach(dataFrame: DG.DataFrame) {
    this.viewSubs = [];
    try {
      // PtmFilter.checkRefDf(this.ptmMap, this.cdrMap, this.refDf);

      console.debug(`MLB: PtmFilter.attach() start, ${((Date.now() - _startInit) / 1000).toString()}`);
      super.attach(dataFrame);
    } catch (err: unknown) {
      console.error(err instanceof Error ? err.message : (err as Object).toString());
    } finally {
      console.debug(`MLB: PtmFilter.attach() end, ${((Date.now() - _startInit) / 1000).toString()}`);
    }
  }

  applyState(state: any) {
    super.applyState(state);

    // OnPropertyChanged() is not called when setting a property in TableView.filters() argument,
    // but applyState() is called with that value
    if ('currentCdr' in state) {
      const value: string = state['currentCdr'];
      const key: string = this.predictedCdrKeys.find((v) => v.toUpperCase() == value.toUpperCase());
      console.debug(`MLB: PtmFilter.applyState(${MlbEvents.CdrChanged}) value ="${value}" -> key="${key}".`);
      this.currentCdr = key ?? 'None';
    }
    if ('predictedPtm' in state)
      this.predictedPtmDf = DG.DataFrame.fromCsv(state['predictedPtm']);

    if ('observedPtm' in state)
      this.observedPtmDf = DG.DataFrame.fromCsv(state['observedPtm']);

    // Calculating indexes is possible only after required data applied
    const predictedTempDf = this.predictedPtmDf.clone(null, ['v_id']);
    predictedTempDf.columns.addNewInt('index').init((i) => i);
    this.predictedIndexes = this.dataFrame.clone(null, ['v id'])
      .join(predictedTempDf, ['v id'], ['v_id'], [], ['index'], 'left', false)
      .getCol('index').getRawData() as Int32Array;

    const observedTempDf = this.observedPtmDf.clone(null, ['v_id']);
    observedTempDf.columns.addNewInt('index').init((i) => i);
    this.observedIndexes = this.dataFrame.clone(null, ['v id'])
      .join(observedTempDf, ['v id'], ['v_id'], [], ['index'], 'left', false)
      .getCol('index').getRawData() as Int32Array;

    this.render();
  }

  detach() {
    super.detach();
    console.log('MLB: PtmFilter.detach() filter detached');

    this.viewSubs.forEach((u) => { u.unsubscribe(); });
    this.viewSubs = null;
  }

  applyFilter() {
    const t1: number = Date.now();
    try {
      console.debug('MLB: PtmFilter.applyFilter() start');

      const getStateFor = (index: number) => {
        for (const checkedPredictedPtm of this.predictedPtmInputValue) {
          const ptmSymbol = this.predictedPtmMap[`${this.chainType} ${checkedPredictedPtm}`];
          const refStr = this.predictedPtmDf.get(ptmSymbol, this.predictedIndexes[index]);
          if (typeof refStr === 'undefined' || refStr === '')
            return false;

          const refVal = JSON.parse(refStr);
          // const refMax: number = refVal['max'];

          const queryCdrName: string = `${this.currentCdr.toLowerCase()}_CDR${this.chainType}_ranges`;
          if (!(queryCdrName in this.predictedCdrMap)) {
            console.warn(`Missing CDR '${queryCdrName}' in CDR map of predicted PTMs.`);
            return false;
          }

          const queryCdrSymbol = this.predictedCdrMap[queryCdrName];
          if (!(queryCdrSymbol in refVal))
            return false;

          const currentProbability = refVal[queryCdrSymbol];
          if (
            typeof currentProbability === 'undefined' ||
            currentProbability < this.predictedProbabilityMin || currentProbability > this.predictedProbabilityMax
          ) return false;
        }

        for (const checkedObservedPtm of this.observedPtmInputValue) {
          const ptmSymbol = this.observedPtmMap[`${this.chainType} ${checkedObservedPtm}`];
          const refStr = this.observedPtmDf.get(ptmSymbol, this.observedIndexes[index]);
          if (!refStr) return false;


          //const refVal = JSON.parse(refStr);
          const refVal = Object.assign({},
            ...refStr.split(';').map((p) => {
              const pp = p.split(':');
              return ({[pp[0]]: pp[1]});
            }));
          // const refMax: number = refVal['o_max'];

          const queryCdrName: string = `${this.currentCdr.toLowerCase()}_CDR${this.chainType}_ranges`;
          if (!(queryCdrName in this.observedCdrMap)) {
            console.warn(`Missing CDR '${queryCdrName}' in CDR map of observed PTMs.`);
            return false;
          }

          const queryCdrSymbol = this.observedCdrMap[queryCdrName];
          if (!(queryCdrSymbol in refVal))
            return false;

          const currentProbability = refVal[queryCdrSymbol];
          if (
            typeof currentProbability === 'undefined' ||
            currentProbability < this.observedProbabilityMin || this.observedProbabilityMax < currentProbability
          ) return false;
        }

        return true;
      };

      const filter = this.dataFrame.filter;
      const rowCount = this.dataFrame.rowCount;
      for (let i = 0; i < rowCount; i++)
        filter.set(i, filter.get(i) && getStateFor(i), false);
      filter.fireChanged();

      console.debug('MLB: PtmFilter.applyFilter() end');
    } catch (err: unknown) {
      console.error(err instanceof Error ? err.message : (err as Object).toString());
    } finally {
      const t2 = Date.now();
      console.debug('MLB: PtmFilter.applyFilter() ET, ' + `${((t2 - t1) / 1000).toString()} sec.`);
    }
  }

  render() {
    try {
      console.debug(`MLB: PtmFilter.render() start, ${((Date.now() - _startInit) / 1000).toString()}`);

      $(this.root).empty();

      const chainTypeInput = ui.choiceInput('Chain type', 'L', ['L', 'H'], () => {
        this.chainType = chainTypeInput.stringValue as ChainTypeType;
        this.dataFrame.rows.requestFilter();
      });

      const predictedCdrInput = ui.choiceInput('CDR', this.currentCdr, this.predictedCdrKeys, () => {
        this.currentCdr = predictedCdrInput.stringValue === 'None' ?
          'max' : this.predictedCdrMap[this.normalizedCDRMap[`${this.chainType} ${predictedCdrInput.stringValue}`]];
        this.dataFrame.rows.requestFilter();
      });

      const [predictedPtmInputLabel, predictedProbabilityHost, predictedPtmInput] = this.buildPredictedPtmInput();
      const [observedPtmInputLabel, observedProbabilityHost, observedPtmInput] = this.buildObservedPtmInput();

      // const ptmInputGrid = this.buildPredictedPtmInputGrid();

      this.root = ui.divV([
        chainTypeInput.root,
        predictedCdrInput.root,
        predictedPtmInputLabel, predictedProbabilityHost, predictedPtmInput.root,
        ui.element('hr'),
        observedPtmInputLabel, observedProbabilityHost, observedPtmInput.root,
        ui.element('hr'),
        // predictedPtmInputLabel, ptmInputGrid.root,
        ui.element('hr'),
      ]);
      this.root.style.margin = '10px';
    } catch (err: unknown) {
      console.error(err instanceof Error ? err.message : (err as Object).toString());
    } finally {
      console.debug(`MLB: PtmFilter.render() end, ${((Date.now() - _startInit) / 1000).toString()}`);
    }
  }

  static adjustPtmMultiChoiceInput(root: HTMLElement) {
    const cash = $(root);

    root.style.marginLeft = '10px';

    cash.find('div > div.ui-input-multi-choice-checks > div > input.ui-input-editor[type="checkbox"]')
      .each((_i, element) => {
        // element.style.margin = '0';
      });
    cash.find('div > div.ui-input-multi-choice-checks > div')
      .each((_i, element) => {
        ui.tooltip.bind(element, `${element.lastChild.textContent}`);
      });
    cash.find('div > div.ui-input-multi-choice-checks > label')
      .each((_i, element) => {
        element.style.marginLeft = '10px';
        ui.tooltip.bind(element, `${element.lastChild.textContent}`);
      });
    // adjust PTM choice inputs width
    cash.find('div > div.ui-input-multi-choice-checks')
      .each((_, el) => {
        // el.style.width = '100%';
        // el.style.width = 'max-content';
        el.style.maxHeight = 'none'; // removes scroller and shows all items of PTM
        el.style.setProperty('padding', '0', 'important');
        el.style.setProperty('overflow-y', 'hidden', 'important'); // disables vertical scroll
        el.style.setProperty('overflow-x', 'hidden', 'important'); // disables horizontal scroll
      });
    cash.find('div > div.ui-input-multi-choice-checks > div > label')
      .each((_, el) => {
        el.style.whiteSpace = 'nowrap';
        el.style.marginLeft = '8px';
        el.style.marginTop = '2px';
      });

    // remove input label
    cash.children().first().remove(); // Remove label of multiChoiceInput 'PTM'
  }

  static adjustPtmProbabilityInput(ptmProbabilityInput: DG.RangeSlider) {
    const probabilityInputCash = $(ptmProbabilityInput.root);
    probabilityInputCash.find('div > svg')
      .each((_, el) => {
        el.style.width = '100%';
        el.parentElement.style.width = '100%';
        el.style.height = '17px';
        el.style.marginTop = '3px';
      });
  }

  static buildProbabilityInput(
    probabilityText: string, cMin: number, cMax: number,
    onChanged: (min: number, max: number, header: HTMLElement) => void
  ): [HTMLElement, Subscription] {
    const header = ui.divText(probabilityText);
    header.style.marginTop = '7px';
    header.style.marginBottom = '5px';
    const input = ui.rangeSlider(0, 1, cMin, cMax, false, 'thin_barbell');

    PtmFilter.adjustPtmProbabilityInput(input);

    const rsMin = ui.inlineText(['0']);
    rsMin.style.marginRight = '5px';
    const rsMax = ui.inlineText(['1']);
    rsMax.style.marginLeft = '3px';
    const host = ui.divV([
      header,
      ui.divH([rsMin, input.root, rsMax]),
    ], {style: {marginLeft: '10px'}});

    const subs: Subscription = input.onValuesChanged.subscribe(() => {
      onChanged(parseFloat(input.min.toFixed(3)), parseFloat(input.max.toFixed(3)), header);
    });

    return [host, subs];
  }

  private buildPredictedPtmInput(): [HTMLLabelElement, HTMLElement, DG.InputBase<string[]>] {
    const predictedPtmInputLabel = ui.label('PTM Predicted');
    predictedPtmInputLabel.style.marginTop = '5px';

    const predictedPtmInput = ui.multiChoiceInput('PTM', this.predictedPtmInputValue, this.predictedPtmKeys, () => {
      this.predictedPtmInputValue = predictedPtmInput.value;
      this.dataFrame.rows.requestFilter();
    });

    PtmFilter.adjustPtmMultiChoiceInput(predictedPtmInput.root);

    const [predictedProbabilityHost, predictedProbabilitySub] = PtmFilter.buildProbabilityInput(
      this.getPredictedPtmProbabilityText(), this.predictedProbabilityMin, this.predictedProbabilityMax,
      (min, max, header) => {
        this.predictedProbabilityMin = min;
        this.predictedProbabilityMax = max;
        header.textContent = this.getPredictedPtmProbabilityText();
        this.dataFrame.rows.requestFilter();
      });
    this.viewSubs.push(predictedProbabilitySub);

    return [predictedPtmInputLabel, predictedProbabilityHost, predictedPtmInput];
  }

  private buildObservedPtmInput(): [HTMLLabelElement, HTMLElement, DG.InputBase<string[]>] {
    const observedPtmInputLabel = ui.label('PTM Observed');
    observedPtmInputLabel.style.marginTop = '5xpx';

    const observedPtmInput = ui.multiChoiceInput('PTM', this.observedPtmInputValue, this.observedPtmKeys, () => {
      this.observedPtmInputValue = observedPtmInput.value;
      this.dataFrame.rows.requestFilter();
    });

    PtmFilter.adjustPtmMultiChoiceInput(observedPtmInput.root);

    const [observedProbabilityHost, observedProbabilitySub] = PtmFilter.buildProbabilityInput(
      this.getObservedPtmProbabilityText(), this.predictedProbabilityMin, this.predictedProbabilityMax,
      (min, max, header) => {
        this.observedProbabilityMin = min;
        this.observedProbabilityMax = max;
        header.textContent = this.getObservedPtmProbabilityText();
        this.dataFrame.rows.requestFilter();
      });
    this.viewSubs.push(observedProbabilitySub);

    return [observedPtmInputLabel, observedProbabilityHost, observedPtmInput];
  }

  private buildPredictedPtmInputGrid() {
    const ptmInputRowHeight = 20;
    const ptmInputDf: DG.DataFrame = DG.DataFrame.fromObjects(this.predictedPtmKeys.map((ptmKey) => {
      return {
        predicted: false, observed: false, name: ptmKey,
        lightSymbol: this.predictedPtmMap[`L ${ptmKey}`], heavySymbol: this.predictedPtmMap[`H ${ptmKey}`]
      };
    }));

    const ptmInputGrid: DG.Grid = ptmInputDf.plot.grid({
      allowEdit: true,
      allowRowSelection: false,
      allowRowResizing: false,
      allowRowReordering: false,
      allowColSelection: false,
      allowColResizing: false,
      allowColReordering: false,
      allowColHeaderResizing: false,
      allowBlockSelection: false,
      showRowHeader: false,
      showColumnGridlines: false,
      showRowGridlines: false,
      // menus
      allowColumnMenu: false,
      showDefaultPopupMenu: false,
      showContextMenu: false,
      topLevelDefaultMenu: false,
      // current
      showCurrentRowIndicator: false,
      showCurrentCellOutline: false,
      // mouse over
      showMouseOverRowIndicator: false,
      showMouseOverRowStripe: false,

      marginLeft: 5,
      rowHeight: ptmInputRowHeight,
    });
    // @formatter:off
    ((col) => { col.name = 'P'; col.width = 20; })(ptmInputGrid.columns.byIndex(1));
    ((col) => { col.name = 'O'; col.width = 20; })(ptmInputGrid.columns.byIndex(2));
    ((col) => { col.name = 'Modification'; col.editable = false; })(ptmInputGrid.columns.byIndex(3));
    ((col) => { col.visible = false; })(ptmInputGrid.columns.byIndex((4))); // lightSymbol
    ((col) => { col.visible = false; })(ptmInputGrid.columns.byIndex((5))); // heavySymbol
    // @formatter:on
    ptmInputGrid.autoSize(450, 11 + ptmInputRowHeight * (this.predictedPtmKeys.length + 1), 25, 25, true);

    this.viewSubs.push(ptmInputGrid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader && ['P', 'O'].includes(cell.gridColumn.name)) {
        let txt = '';
        switch (cell.gridColumn.name) {
        case 'P':
          txt = 'Predicted';
          break;
        case 'O':
          txt = 'Observed';
          break;
        }
        ui.tooltip.show(txt, x + 8, y + 8);
      }
      return true;
    }));
    this.viewSubs.push(ptmInputGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => {
      if (args.cell.isColHeader) {
        if (args.cell.gridColumn.visible) {
          args.g.font = args.g.font.replace('bold', '');
          args.g.textAlign = ['P', 'O'].includes(args.cell.gridColumn.name) ? 'center' : 'left';
          args.g.fillStyle = '#4b4b4a';
          const textSize = args.g.measureText(args.cell.gridColumn.name);

          let txtX: number;
          switch (args.g.textAlign) {
          case 'center':
            txtX = 4 + args.bounds.x + (args.bounds.width - textSize.width) / 2;
            break;
          case 'left':
            txtX = 8 + args.bounds.x;
            break;
          }
          const txtY: number = args.bounds.y +
            (textSize.fontBoundingBoxAscent /* + textSize.fontBoundingBoxDescent /**/);
          args.g.fillText(args.cell.gridColumn.name, txtX, txtY, args.bounds.width);
        }
        args.preventDefault(); // this is required to prevent drawing headers of hidden columns
      }
    }));
    ptmInputDf.onDataChanged.subscribe((value) => {
      this.predictedPtmInputValue = wu.count().take(ptmInputDf.rowCount)
        .map((rowI) => [ptmInputDf.get('predicted', rowI), ptmInputDf.get('name', rowI)])
        .filter((v) => v[0]).map((v) => v[1]).toArray();
      this.observedPtmInputValue = wu.count().take(ptmInputDf.rowCount)
        .map((rowI) => [ptmInputDf.get('observed', rowI), ptmInputDf.get('name', rowI)])
        .filter((v) => v[0]).map((v) => v[1]).toArray();
      this.dataFrame.rows.requestFilter();
    });

    return ptmInputGrid;
  }

  private getPredictedPtmProbabilityText(): string {
    return `Probability: [${this.predictedProbabilityMin}; ${this.predictedProbabilityMax}]`;
  }

  private getObservedPtmProbabilityText(): string {
    return `Probability: [${this.observedProbabilityMin}; ${this.observedProbabilityMax}]`;
  }
}
