import * as DG from 'datagrok-api/dg';
import {DEFAULT_FORMATS} from '../const';
import {GROUP_TYPE, PHOSPHATE_SYMBOL, UNKNOWN_SYMBOL} from './const';
import {CodesInfo} from '../data-loading-utils/types';
import {codesToHelmDictionary} from '../data-loading-utils/json-loader';

const HELM_WRAPPER = {
  LEFT: 'RNA1{',
  RIGHT: '}$$$$',
};

export class FormatConverter {
  constructor(private readonly sequence: string, private readonly sourceFormat: string) { };

  convertTo(targetFormat: string): string {
    const formats = Object.keys(codesToHelmDictionary);

    if (this.sourceFormat === DEFAULT_FORMATS.HELM && formats.includes(targetFormat))
      return helmToFormat(this.sequence, targetFormat);
    else if (formats.includes(this.sourceFormat) && targetFormat === DEFAULT_FORMATS.HELM)
      return formatToHelm(this.sequence, this.sourceFormat);
    else if ([this.sourceFormat, targetFormat].every((el) => formats.includes(el))) {
      const helm = formatToHelm(this.sequence, this.sourceFormat);
      return helmToFormat(helm, targetFormat);
    }
    else {
      throw new Error (`ST: unsupported translation direction ${this.sourceFormat} -> ${targetFormat}`);
    }
  }
}

function getRegExpPattern(arr: string[]): string {
  const negativeLookBehind = '(?<!\\([^()]*)'; // not '(' followed by non-parenths
  const negativeLookAhead = '(?![^()]*\\))';  // not ')' preceded by non-parenths
  const escaped = arr.map((key) => key.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'))
    .map((key) => {
    if (!key.includes('(') && !key.includes(')'))
      return `${negativeLookBehind}${key}${negativeLookAhead}`;
    return key;
  });
  const result =  escaped.join('|');
  return result;
}

function sortCallback(a: string, b: string) {return b.length - a.length};

function getHelmToCodeDict(infoObj: CodesInfo) {
  const result: {[key: string]: string | string[]} = {};
  Object.values(infoObj).forEach((obj: {[code: string]: string}) => {
    Object.entries(obj).forEach(([code, helm]) => {
      const key = helm.replace(/\)p/g, ')').replace(/\]p/g, ']');
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

function helmToFormat(helmSequence: string, targetFormat: string): string {
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
  }).replace(/\?+/g, UNKNOWN_SYMBOL).replace(/p\.|\./g, '');
  result = result.replace(/<empty>/g, '');
  // remove double slash in LCMS codes
  result = result.replace(/\/\//g, '/');
  return result;
}

function formatToHelm(sequence: string, sourceFormat: string): string {
  const codesInfoObject = codesToHelmDictionary[sourceFormat] as CodesInfo;
  const dict = Object.assign({}, ...Object.values(codesInfoObject)) as {[code: string]: string};

  const formatCodes = Object.keys(dict).sort(sortCallback);
  const formatRegExp = new RegExp(getRegExpPattern(formatCodes) + '|\\([^()]*\\)|.', 'g'); // the added group before '|.' is to avoid mismatch inside parenths

  const phosphateHELMCodes = Array.from(
    new Set(Object.values(codesInfoObject[GROUP_TYPE.LINKAGE]))
  ).sort(sortCallback);
  const phosphateHELMPattern = getRegExpPattern(phosphateHELMCodes);
  const phosphateRegExp = new RegExp(`${PHOSPHATE_SYMBOL}\.(${phosphateHELMPattern})`, 'g');

  let helm = sequence.replace(formatRegExp, (match) => {
    const result = formatCodes.includes(match) ?  dict[match] + '.' : '?';
    return result;
  });
  helm = helm.replace(/\?+/g, `${UNKNOWN_SYMBOL}.`);
  helm = helm.slice(0, -1); // strip last dot
  if (helm[helm.length - 1] === PHOSPHATE_SYMBOL)
    helm = helm.slice(0, -1);
  helm = helm.replace(phosphateRegExp, (match, group) => group);
  helm = helm.replace(/<empty>/g, '');
  return `${HELM_WRAPPER.LEFT + helm + HELM_WRAPPER.RIGHT}`;
}
