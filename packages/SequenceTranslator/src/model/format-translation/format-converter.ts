import * as DG from 'datagrok-api/dg';
import {DELIMITER} from '../const';
import {sortByReverseLength} from '../helpers';
import {MonomerLibWrapper} from '../monomer-lib/lib-wrapper';
import {SYNTHESIZERS as FORMAT} from '../const';
import {GROUP_TYPE, EDGES, CENTER, PHOSPHORUS} from './const';
import {KeyToValue, CodesInfo} from '../data-loading-utils/types';
import {formatDictionary, codesToHelmDictionary} from '../data-loading-utils/json-loader';

const HELM_WRAPPER = {
  LEFT: 'RNA1{',
  RIGHT: '}$$$$',
};

// todo: remove strange legacy logic with magic numbers
export class FormatConverter {
  constructor(private readonly sequence: string, private readonly sourceFormat: FORMAT) { };

  convertTo(targetFormat: FORMAT): string {
    if (this.sourceFormat === FORMAT.HELM)
      return this.helmToGcrs();
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

  private simpleConversion(codeMapping: KeyToValue) {
    const regexp = new RegExp(getRegExpPattern(sortByReverseLength(Object.keys(codeMapping))), 'g');
    return this.sequence.replace(regexp, (code) => codeMapping[code]);
  }

  private bioSpringToGcrs(codeMapping: KeyToValue): string {
    let count: number = -1;
    const regexp = new RegExp(getRegExpPattern(Object.keys(codeMapping)), 'g');
    return this.sequence.replace(regexp, (x: string) => {
        count++;
        return (count == 4) ? codeMapping[x].slice(0, -3) + 'ps' : (count == 14) ? codeMapping[x].slice(0, -2) + 'nps' : codeMapping[x];
      });
  }

  private gcrsToHelm(): string {
    const codesInfoObject = codesToHelmDictionary[FORMAT.GCRS] as CodesInfo;
    const dict = Object.assign({}, ...Object.values(codesInfoObject)) as {[code: string]: string};

    const gcrsCodes = Object.keys(dict).sort(sortCallback);
    const gcrsRegExp = new RegExp(getRegExpPattern(gcrsCodes), 'g');

    const phosphateHELMCodes = Array.from(
      new Set(Object.values(codesInfoObject[GROUP_TYPE.LINKAGE]))
    ).sort(sortCallback);
    const phosphateHELMPattern = getRegExpPattern(phosphateHELMCodes);
    const phosphateRegExp = new RegExp(`${PHOSPHORUS}\.(${phosphateHELMPattern})`, 'g');

    let helm = this.sequence.replace(gcrsRegExp, (match) => dict[match] + '.');
    helm = helm.slice(0, -1); // strip last dot
    if (helm[helm.length - 1] === PHOSPHORUS) 
      helm = helm.slice(0, -1);
    helm = helm.replace(phosphateRegExp, (match, group) => group);
    return `${HELM_WRAPPER.LEFT + helm + HELM_WRAPPER.RIGHT}`;
  }

  private helmToGcrs(): string {
    const codesInfoObject = codesToHelmDictionary[FORMAT.GCRS] as CodesInfo;
    const dict = getHelmToCodeDict(codesInfoObject);
    const wrapperRegExp = new RegExp(getRegExpPattern(Object.values(HELM_WRAPPER)), 'g')
    let gcrs = this.sequence.replace(wrapperRegExp, '');

    const helmCodes = Object.keys(dict)
      .sort(sortCallback);
    const helmRegExp = new RegExp(getRegExpPattern(helmCodes), 'g');
    gcrs = gcrs.replace(helmRegExp, (match) => dict[match]);
    gcrs = gcrs.replace(/p\.|\./g, '');
    return gcrs;
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
    const regexp = new RegExp(getRegExpPattern(Object.keys(edgeCodeMapping)), 'g')
    return this.sequence.replace(regexp, (x: string) => {
      count++;
      return (count > 4 && count < 15) ? centerCodeMapping[x] : edgeCodeMapping[x];
    });
  }

  private nucleotidesToGCRS(edgeCodeMapping: KeyToValue, centerCodeMapping: KeyToValue): string {
    let count: number = -1;
    const regexp = new RegExp(getRegExpPattern(Object.keys(edgeCodeMapping)), 'g');
    return this.sequence.replace(regexp, (x: string) => {
      count++;
      if (count < 5) return (count == 4) ? edgeCodeMapping[x].slice(0, -3) + 'ps' : edgeCodeMapping[x];
      if (count < 15) return (count == 14) ? centerCodeMapping[x].slice(0, -2) + 'nps' : centerCodeMapping[x];
      return edgeCodeMapping[x];
    });
  }
}

function getRegExpPattern(arr: string[]): string {
  const escaped = arr.map((key) => key.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'));
  return escaped.join('|');
}

function sortCallback(a: string, b: string) {return b.length - a.length};

function getHelmToCodeDict(infoObj: CodesInfo) {
  const result: {[key: string]: string | string[]} = {};
  Object.values(infoObj).forEach((obj: {[code: string]: string}) => {
    Object.entries(obj).forEach(([code, helm]) => {
      const key = helm.replace(/\)p/g, ')');
      if (result[key] === undefined) {
        result[key] = [code];
      } else {
        (result[key] as string[]).push(code);
      }
    })
  });
  Object.entries(result).forEach(([key, value]) => {
    const sorted = (value as string[]).sort(sortCallback);
    result[key] = sorted[0] as string;
  })
  return result as {[key: string]: string};
}
