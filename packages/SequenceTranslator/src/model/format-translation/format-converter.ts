import * as DG from 'datagrok-api/dg';
import {sortByReverseLength} from '../helpers';
import {FORMAT} from '../const';
import {GROUP_TYPE, EDGES, CENTER, PHOSPHATE} from './const';
import {KeyToValue, CodesInfo} from '../data-loading-utils/types';
import {formatDictionary, codesToHelmDictionary} from '../data-loading-utils/json-loader';

const HELM_WRAPPER = {
  LEFT: 'RNA1{',
  RIGHT: '}$$$$',
};

export class FormatConverter {
  constructor(private readonly sequence: string, private readonly sourceFormat: FORMAT) { };

  convertTo(targetFormat: FORMAT): string {
    const formatsConvertibleToHelm = Object.keys(codesToHelmDictionary);

    if (this.sourceFormat === FORMAT.HELM && formatsConvertibleToHelm.includes(targetFormat))
      return helmToFormat(this.sequence, targetFormat);

    if (formatsConvertibleToHelm.includes(this.sourceFormat) && targetFormat === FORMAT.HELM)
      return formatToHelm(this.sequence, this.sourceFormat);

    if ([this.sourceFormat, targetFormat].every((el) => formatsConvertibleToHelm.includes(el))) {
      const helm = formatToHelm(this.sequence, this.sourceFormat);
      return helmToFormat(helm, targetFormat);
    }

    const codeMapping = formatDictionary[this.sourceFormat][targetFormat];

    if (codeMapping === undefined && targetFormat !== FORMAT.HELM) {
      throw new Error (`ST: unsupported translation direction ${this.sourceFormat} -> ${targetFormat}`);
    }

    if (this.sourceFormat === FORMAT.NUCLEOTIDES) {
      const edgeCodeMapping = codeMapping[EDGES] as KeyToValue;
      const centerCodeMapping = codeMapping[CENTER] as KeyToValue;
      if (targetFormat === FORMAT.BIOSPRING) {
        return this.nucleotidesToBioSpring(edgeCodeMapping, centerCodeMapping);
      } else { // target === GCRS
        return this.nucleotidesToGCRS(edgeCodeMapping, centerCodeMapping);
      }
    } else
      return this.gcrsToNucleotides(codeMapping as KeyToValue);
  }

  private gcrsToNucleotides(codeMapping: KeyToValue) {
    const regexp = new RegExp(getRegExpPattern(sortByReverseLength(Object.keys(codeMapping))), 'g');
    return this.sequence.replace(regexp, (code) => codeMapping[code]);
  }

  // todo: remove strange legacy logic with magic numbers
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
      if (count < 5) return (count === 4) ? edgeCodeMapping[x].slice(0, -3) + 'ps' : edgeCodeMapping[x];
      if (count < 15) return (count === 14) ? centerCodeMapping[x].slice(0, -2) + 'nps' : centerCodeMapping[x];
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

function helmToFormat(helmSequence: string, targetFormat: FORMAT): string {
  const codesInfoObject = codesToHelmDictionary[targetFormat] as CodesInfo;
  const dict = getHelmToCodeDict(codesInfoObject);
  const wrapperRegExp = new RegExp(getRegExpPattern(Object.values(HELM_WRAPPER)), 'g')
  let result = helmSequence.replace(wrapperRegExp, '');

  const helmCodes = Object.keys(dict)
    .sort(sortCallback);
  const helmRegExp = new RegExp(getRegExpPattern(helmCodes) + '|.', 'g');
  result = result.replace(helmRegExp, (match) => {
    return helmCodes.includes(match) ? dict[match] :
      (match === 'p' || match === '.') ? match : '?';
  }).replace(/\?+/g, '<?>').replace(/p\.|\./g, '');
  return result;
}

function formatToHelm(sequence: string, sourceFormat: FORMAT): string {
  const codesInfoObject = codesToHelmDictionary[sourceFormat] as CodesInfo;
  const dict = Object.assign({}, ...Object.values(codesInfoObject)) as {[code: string]: string};

  const formatCodes = Object.keys(dict).sort(sortCallback);
  const formatRegExp = new RegExp(getRegExpPattern(formatCodes), 'g');

  const phosphateHELMCodes = Array.from(
    new Set(Object.values(codesInfoObject[GROUP_TYPE.LINKAGE]))
  ).sort(sortCallback);
  const phosphateHELMPattern = getRegExpPattern(phosphateHELMCodes);
  const phosphateRegExp = new RegExp(`${PHOSPHATE}\.(${phosphateHELMPattern})`, 'g');

  let helm = sequence.replace(formatRegExp, (match) => dict[match] + '.');
  helm = helm.slice(0, -1); // strip last dot
  if (helm[helm.length - 1] === PHOSPHATE) 
    helm = helm.slice(0, -1);
  helm = helm.replace(phosphateRegExp, (match, group) => group);
  return `${HELM_WRAPPER.LEFT + helm + HELM_WRAPPER.RIGHT}`;
}
