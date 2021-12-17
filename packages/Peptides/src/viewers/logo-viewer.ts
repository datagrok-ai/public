import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import * as logojs from 'logojs-react';
import {splitAlignedPeptides} from '../utils/split-aligned';
import {ChemPalette} from '../utils/chem-palette';

/**
 * Logo viewer.
 *
 * @export
 * @class Logo
 * @extends {DG.JsViewer}
 */
export class Logo extends DG.JsViewer {
  initialized: boolean;
  option: any;
  colSemType: string;
  splitted: DG.DataFrame | null;
  ppm: Array<Array<number>>;
  reactHost: HTMLDivElement | null;
  PROT_NUMS: { [id: string]: number };
  LET_COLORS: Array<any>;
  target: DG.DataFrame | undefined | null;

  /**
   * Creates an instance of Logo.
   * 
   * @memberof Logo
   */
  constructor() {
    super();
    this.initialized = false;
    this.colSemType = this.string('colSemType', 'alignedSequence');

    this.splitted = null;
    this.ppm = [];
    this.reactHost = null;
    this.target = null;
    this.PROT_NUMS = {
      'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
      'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'U': 18, 'V': 19, 'W': 20, 'Y': 21, 'Z': 22,
    };
    //TODO: use chem palette
    this.LET_COLORS = [
      {color: 'rgb(44,160,44)', regex: 'A'},
      {color: 'rgb(44,160,44)', regex: 'B'},
      {color: 'rgb(188,189,34)', regex: 'C'},
      {color: 'rgb(31,119,180)', regex: 'D'},
      {color: 'rgb(30,110,96)', regex: 'E'},
      {color: 'rgb(24,110,79)', regex: 'F'},
      {color: 'rgb(214,39,40)', regex: 'G'},
      {color: 'rgb(158,218,229)', regex: 'H'},
      {color: 'rgb(122, 102, 189)', regex: 'I'},
      {color: 'rgb(108, 218, 229)', regex: 'K'},
      {color: 'rgb(30,110,96)', regex: 'L'},
      {color: 'rgb(141, 124, 217)', regex: 'M'},
      {color: 'rgb(235,137,70)', regex: 'N'},
      {color: 'rgb(255,152,150)', regex: 'P'},
      {color: 'rgb(205, 111, 71)', regex: 'Q'},
      {color: 'rgb(23,190,207)', regex: 'R'},
      {color: 'rgb(255,187,120)', regex: 'S'},
      {color: 'rgb(245,167,100)', regex: 'T'},
      {color: 'rgb(188,189,34)', regex: 'U'},
      {color: 'rgb(23,190,207', regex: 'V'},
      {color: 'rgb(182, 223, 138)', regex: 'W'},
      {color: 'rgb(152,223,138)', regex: 'Y'},
      {color: 'rgb(205, 111, 71)', regex: 'Z'},
    ];
  }

  /**
   * Initializer function.
   *
   * @memberof Logo
   */
  init() {
    this.initialized = true;
    console.log('INIT');
    this.target = this.dataFrame;
    [this.splitted] = splitAlignedPeptides(this.dataFrame!.columns.bySemType(this.colSemType));
    this.root.style.width = 'auto';
    this.root.style.height = 'auto';
    this.root.style.maxHeight = '200px';
  }

  /**
   * Function to execute when the table is attached.
   *
   * @memberof Logo
   */
  onTableAttached() {
    if (typeof this.dataFrame !== 'undefined') {
      if (!this.initialized) {
        this.init();
      }

      this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_: any) => this.render()));
      this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_: any) => this.render()));
    }

    this.render();
  }

  /**
   * Function that is executed when the viewer is detached.
   *
   * @memberof Logo
   */
  detach() {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /**
   * Function that is executed when the viewer property is changed.
   *
   * @param {DG.Property} property
   * @memberof Logo
   */
  onPropertyChanged(property: DG.Property) {
    super.onPropertyChanged(property);

    this.render();
  }

  /**
   * Function that renders the viewer.
   *
   * @memberof Logo
   */
  async render() {
    const bits = this.dataFrame!.selection;
    let selected = false;
    if (bits.trueCount > 0) {
      selected = true;
      this.target = this.dataFrame!
        .groupBy([this.dataFrame!.columns.bySemType(this.colSemType).name])
        .whereRowMask(this.dataFrame!.selection)
        .aggregate();
    }
    if (selected) {
      [this.splitted] = splitAlignedPeptides(this.target!.columns.bySemType(this.colSemType));
    } else [this.splitted] = splitAlignedPeptides(this.dataFrame!.columns.bySemType(this.colSemType));
    $(this.root).empty();

    if (typeof this.dataFrame !== 'undefined') {
      this.findLogo();
    }
  }

  /**
   * Create logo.
   *
   * @memberof Logo
   */
  async findLogo() {
    this.getInfoFromDf();
    logojs.embedProteinLogo(this.root, {alphabet: this.LET_COLORS, ppm: this.ppm});
  }

  /**
   * Retrieves information for building logo from the dataframe.
   *
   * @memberof Logo
   */
  getInfoFromDf() {
    let index: number = 0;
    this.ppm = [];

    for (const col of this.splitted!.columns) {
      const size = col.length;
      this.ppm.push(new Array(22).fill(0));
      for (let i = 0; i < col.length; i++) {
        const c = col.get(i);
        if (c != '-') {
          if (c[1] == '(') {
            this.ppm[index][this.PROT_NUMS[c.substr(0, 1).toUpperCase()]] += 1 / size;
          } else if (c.substr(0, 3) in ChemPalette.AAFullNames && (c.length == 3 || c.at(3) == '(')) {
            this.ppm[index][this.PROT_NUMS[ChemPalette.AAFullNames[c.substr(0, 3)]]] += 1 / size;
          } else if (c.at(0)?.toLowerCase() == c.at(0) && c.substr(1, 3) in ChemPalette.AAFullNames &&
            (c.length == 4 || c.at(4) == '(')
          ) {
            this.ppm[index][this.PROT_NUMS[ChemPalette.AAFullNames[c.substr(1, 3)]]] += 1 / size;
          } else {
            this.ppm[index][this.PROT_NUMS[c]] += 1 / size;
          }
        }
      }
      index++;
    }
  }
}
