import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import $ from 'cash-dom';
import * as logojs from 'logojs-react';

export class Logo extends DG.JsViewer {
  initialized: boolean;
  option: any;
  colSemType: string;
  splitted: DG.DataFrame | null;
  ppm: Array<Array<number>>;
  df: DG.DataFrame;

  LET_COLORS = [
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

  PROT_NUMS: { [name: string]: number } = {
    'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,
    'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'U': 18, 'V': 19, 'W': 20, 'Y': 21, 'Z': 22,
  };

  constructor(splittedTable: DG.DataFrame, mlbTable: DG.DataFrame) {
    super();
    this.initialized = false;
    this.colSemType = this.string('colSemType', 'alignedSequence');

    this.splitted = null;
    this.ppm = [];

    this.initialized = true;
    this.df = mlbTable;
    this.splitted = splittedTable;

    this.root.style.width = 'auto';
    this.root.style.height = 'auto';
    this.root.style.maxHeight = '200px';
  }


  render(newAligned: DG.DataFrame) {
    this.splitted = newAligned;
    $(this.root).empty();
    this.findLogo();
  }

  findLogo() {
    this.getInfoFromDf();
    logojs.embedProteinLogo(this.root, {alphabet: this.LET_COLORS, ppm: this.ppm});
  }

  getInfoFromDf() {
    let index: number = 0;
    this.ppm = [];

    for (const col of this.splitted!.columns) {
      const size = col.length;
      this.ppm.push(new Array(22).fill(0));
      for (let i = 0; i < col.length; i++) {
        const c = col.get(i);
        if (c != '-') {
          if (c[1] == '(')
            this.ppm[index][this.PROT_NUMS[c.substr(0, 1).toUpperCase()]] += 1 / size;
          else if (c.substr(0, 3) in bio.Aminoacids.AAFullNames && (c.length == 3 || c.at(3) == '('))
            this.ppm[index][this.PROT_NUMS[bio.Aminoacids.AAFullNames[c.substr(0, 3)]]] += 1 / size;
          else if (c.at(0)?.toLowerCase() == c.at(0) && c.substr(1, 3) in bio.Aminoacids.AAFullNames &&
            (c.length == 4 || c.at(4) == '(')
          )
            this.ppm[index][this.PROT_NUMS[bio.Aminoacids.AAFullNames[c.substr(1, 3)]]] += 1 / size;
          else
            this.ppm[index][this.PROT_NUMS[c]] += 1 / size;
        }
      }
      index++;
    }
  }
}
