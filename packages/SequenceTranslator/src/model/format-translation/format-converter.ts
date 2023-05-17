import * as DG from 'datagrok-api/dg';
import {DELIMITER} from '../const';
import {sortByReverseLength} from '../helpers';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {SYNTHESIZERS as FORMAT} from '../const';
import {KeyToValue} from '../data-loading-utils/types';
import {formatDictionary} from '../data-loading-utils/json-loader';

const EDGES = 'edges';
const CENTER = 'center';

// todo: remove strange legacy logic with magic numbers
export class FormatConverter {
  constructor(private sequence: string, private sourceFormat: FORMAT) { };

  convertTo(targetFormat: FORMAT): string {
    const codeMapping = formatDictionary[this.sourceFormat][targetFormat];
    if (codeMapping === undefined && targetFormat !== FORMAT.HELM) {
      throw new Error (`ST: unsupported translation direction ${this.sourceFormat} -> ${targetFormat}`);
    } else if (this.sourceFormat === FORMAT.BIOSPRING && targetFormat === FORMAT.GCRS)
      return this.bioSpringToGcrs(codeMapping as KeyToValue);
    else if (this.sourceFormat === FORMAT.GCRS && targetFormat === FORMAT.LCMS)
      return this.gcrsToLcms(codeMapping as KeyToValue);
    else if (this.sourceFormat === FORMAT.GCRS && targetFormat === FORMAT.HELM)
      return this.gcrsToHelm();
    else if (this.sourceFormat === FORMAT.NUCLEOTIDES) {
      const edgeCodeMapping = codeMapping[EDGES] as KeyToValue;
      const centerCodeMapping = codeMapping[CENTER] as KeyToValue;
      if (targetFormat === FORMAT.BIOSPRING) {
        return this.nucleotidesToBioSpring(edgeCodeMapping, centerCodeMapping);
      } else { // target === GCRS
        return this.nucleotidesToGCRS(edgeCodeMapping, centerCodeMapping);
      }
    } else
      return this.simpleConversion(codeMapping as KeyToValue);
  }

  private buildRegex(keys: string[]): RegExp {
    const escaped = keys.map(key => key.replace(/[.*+?^${}()|\\]/g, '\\$&'));
    return new RegExp(escaped.join('|'), 'g');
  }

  private simpleConversion(codeMapping: KeyToValue) {
    const regex = this.buildRegex(sortByReverseLength(Object.keys(codeMapping)));
    return this.sequence.replace(regex, (code) => codeMapping[code]);
  }

  private bioSpringToGcrs(codeMapping: KeyToValue): string {
    let count: number = -1;
    return this.sequence.replace(this.buildRegex(Object.keys(codeMapping)), (x: string) => {
        count++;
        return (count == 4) ? codeMapping[x].slice(0, -3) + 'ps' : (count == 14) ? codeMapping[x].slice(0, -2) + 'nps' : codeMapping[x];
      });
  }

  private gcrsToHelm(): string {
    function getPattern(arr: string[]): string {
      const escaped = arr.map((key) => key.replace(/[.*+?^${}()|\\]/g, '\\$&'));
      return escaped.join('|');
    }
    const lib = MonomerLibWrapper.getInstance();
    const nucleotideCodes = lib.getGcrsCodesWithoutLinkages().sort((a, b) => b.length - a.length);
    const phosphateCodes = lib.getCodesByFormat(FORMAT.GCRS).filter((el) => !nucleotideCodes.includes(el)).sort((a, b) => b.length - a.length);
    const nucleotideRegExp = new RegExp(getPattern(nucleotideCodes), 'g');
    const phosphateRegExp = new RegExp(`P.(${getPattern(phosphateCodes)})`, 'g');
    this.sequence = this.sequence.replace(nucleotideRegExp, (match) => {
      let group = match;
      if (group.length > 1) {
        group = `[${group}]`;
      }
      group = group.replace(/[()]/g, '');
      return `D(${group})P.`;
    })
      .replace(phosphateRegExp, (match, group) => `[${group}].`)
      .slice(0, -1);
    if (this.sequence[this.sequence.length - 1] === 'P') 
      this.sequence = this.sequence.slice(0, -1);
    return `DNA1{${this.sequence}}$$$`;
    // gcrs.sort((a, b) => b.length - a.length);
    // nucleotides.sort((a, b) => b.length - a.length);
    // const escaped = nucleotides.map((key) => key.replace(/[.*+?^${}()|\\]/g, '\\$&'));
    // const pattern = escaped.join('|');
    // console.log('pattern:', pattern);
    // const nucleotideCode = new RegExp(pattern, 'g');
    // sequence = sequence.replace(nucleotideCode, (match) => {
    //   console.log('match:', match);
    //   let group = match;
    //   if (group.length > 1) {
    //     group = `[${group}]`;
    //   }
    //   group = group.replace(/[()]/g, '');
    //   return `D(${group})P.`;
    // })
    //   .replace(/P.ps/g, '[ps].')
    //   .slice(0, -1);
    // return `DNA1{${sequence}}$$$`;
  }

  // private helmToGcrs(codeMapping: KeyToValue): string {
  // }

  private gcrsToLcms(codeMapping: KeyToValue): string {
    try {
      const lib = MonomerLibWrapper.getInstance();
      codeMapping[DELIMITER] = DELIMITER;
      const codes = Object.keys(codeMapping)
        .concat(DELIMITER)
        .concat(lib.getModificationGCRSCodes());
      const sortedCodes = sortByReverseLength(codes);
      let i = 0;
      let r1 = '';
      while (i < this.sequence.length) {
        const matchedCode = sortedCodes.find((c) => c == this.sequence.slice(i, i + c.length))!;
        r1 += codeMapping[this.sequence.slice(i, i + matchedCode.length)];
        i += matchedCode.length;
      }
      while (r1.indexOf('//') != -1)
        r1 = r1.replace('//', '/');
      return r1;
    } catch {
      return '<error>';
    }
  }

  private nucleotidesToBioSpring(edgeCodeMapping: KeyToValue, centerCodeMapping: KeyToValue): string {
    let count: number = -1;
    return this.sequence.replace(this.buildRegex(Object.keys(edgeCodeMapping)), (x: string) => {
      count++;
      return (count > 4 && count < 15) ? centerCodeMapping[x] : edgeCodeMapping[x];
    });
  }

  private nucleotidesToGCRS(edgeCodeMapping: KeyToValue, centerCodeMapping: KeyToValue): string {
    let count: number = -1;
    return this.sequence.replace(this.buildRegex(Object.keys(edgeCodeMapping)), (x: string) => {
      count++;
      if (count < 5) return (count == 4) ? edgeCodeMapping[x].slice(0, -3) + 'ps' : edgeCodeMapping[x];
      if (count < 15) return (count == 14) ? centerCodeMapping[x].slice(0, -2) + 'nps' : centerCodeMapping[x];
      return edgeCodeMapping[x];
    });
  }
}
