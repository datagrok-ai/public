import * as DG from 'datagrok-api/dg';
import {DELIMITER} from '../const';
import {sortByReverseLength} from '../helpers';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {SYNTHESIZERS as FORMAT} from '../const';
import {KeyToValue} from '../data-loading-utils/types';
import {formatDictionary, codesToHelmDictionary} from '../data-loading-utils/json-loader';

const EDGES = 'edges';
const CENTER = 'center';
const LINKAGE = 'linkage';
const PHOSPHORUS = 'p';

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
      const escaped = arr.map((key) => key.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'));
      return escaped.join('|');
    }

    function sortCallback(a: string, b: string) {return b.length - a.length};

    const dict = codesToHelmDictionary[FORMAT.GCRS] as {[key: string]: string[]};

    const gcrsCodes = Object.keys(dict).sort(sortCallback);
    const gcrsRegExp = new RegExp(getPattern(gcrsCodes), 'g');

    this.sequence = this.sequence.replace(gcrsRegExp, (match) => dict[match][0] + '.');

    const phosphateCodesSet = new Set<string>();
    gcrsCodes.forEach((code) => {
      if (dict[code].includes(LINKAGE))
        phosphateCodesSet.add(dict[code][0]);
    });
    const phosphateCodes = Array.from(phosphateCodesSet).sort(sortCallback);
    const phosphateCodesPattern = getPattern(phosphateCodes.map((el) => el.replaceAll('.', '')));
    const phosphateRegExp = new RegExp(`${PHOSPHORUS}.(${phosphateCodesPattern})`, 'g');

    this.sequence = this.sequence.slice(0, -1); // strip last dot
    if (this.sequence[this.sequence.length - 1] === PHOSPHORUS) 
      this.sequence = this.sequence.slice(0, -1);

    this.sequence = this.sequence.replace(phosphateRegExp, (match, group) => group);
    return `RNA1{${this.sequence}}$$$$`;
  }

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
