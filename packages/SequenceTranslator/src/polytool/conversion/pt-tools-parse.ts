import {getInnerIdx, getOuterIdx, Linkage} from './pt-misc';
import {cleanupHelmSymbol} from '@datagrok-libraries/bio/src/helm/utils';
import {getMonomerPairs, RuleLink, RuleReaction, Rules} from './pt-rules';

type LinkedPosition = {
  firstIdx: number,
  secondIdx: number,
  ruleIdx: number
}

export function parseSeparator(sequence: string): [Linkage[], string[][]] {
  const mainFragments: string[][] = [];
  const linkages: Linkage[] = [];

  const rxp = /(\(.\d+\))?\{[^\}]*\}/g;
  const seqs: string [] = [];
  seqs.push(sequence.replaceAll(rxp, ''));

  //const l = (rxpRes?.length) ?? -1;

  const matches = sequence.matchAll(rxp);
  //const rxpRes = rxp.exec(sequence);
  for (const m of matches) {
    const str = m![0];
    if (str)
      seqs.push(str);
  }

  let counter = 0;
  for (let i = 0; i < seqs.length; i++) {
    const splMonomers = seqs[i].split('-');
    const monomers: string [] = new Array<string>(splMonomers.length);
    let spmCount: number = 0;
    for (let j = 0; j < splMonomers.length; j++) {
      const monomer = splMonomers[j].replace('{', '').replace('}', '');
      if (monomer !== '') {
        monomers[j] = monomer;
        counter++;
        spmCount++;
      } else {
        linkages.push({fChain: i, sChain: i + 1, fMonomer: counter, sMonomer: counter + 1, fR: 1, sR: 1});
      }
    }
    mainFragments.push(monomers.slice(0, spmCount));
  }

  return [linkages, mainFragments];
}

export function parseHelm(sequence: string): [Linkage[], string[][]] {
  const ea = /(\w+\{.*\})\$(.*)\$(.*)\$(.*)\$/g.exec(sequence)!;
  // const fragmentation = helm.split('$');
  const fragmentation = [ea[1], ea[2], ea[3], ea[4]];

  const rawFragments = fragmentation[0].split('|');
  const rawLinkages = fragmentation[1].split('|');

  const monomers = new Array<Array<string>>(rawFragments.length);
  const linkages: Linkage[] = [];

  //HELM parsing
  for (let i = 0; i < rawFragments.length; i++) {
    const idxStart = rawFragments[i].indexOf('{');
    const idxEnd = rawFragments[i].indexOf('}');

    monomers[i] = rawFragments[i].slice(idxStart + 1, idxEnd).split('.').map((s) => cleanupHelmSymbol(s));
  }

  for (let i = 0; i < rawLinkages.length; i++) {
    if (rawLinkages[i] !== '' && rawLinkages[i] !== 'V2.0') {
      const rawData = rawLinkages[i].split(',');
      const fChainIdx = parseInt(rawData[0].replace('PEPTIDE', '')) - 1;
      const sChainIdx = parseInt(rawData[1].replace('PEPTIDE', '')) - 1;
      const rawDataConnections = rawData[2].split('-');
      const rawDataConnection1 = rawDataConnections[0].split(':');
      const rawDataConnection2 = rawDataConnections[1].split(':');

      linkages.push({
        fChain: fChainIdx,
        sChain: sChainIdx,
        fMonomer: getOuterIdx(parseInt(rawDataConnection1[0]), fChainIdx, monomers),
        sMonomer: getOuterIdx(parseInt(rawDataConnection2[0]), sChainIdx, monomers),
        fR: parseInt(rawDataConnection1[1].replace('R', '')),
        sR: parseInt(rawDataConnection2[1].replace('R', '')),
      });
    }
  }

  return [linkages, monomers];
}

export function fromObjectsToHelm(linkages: Linkage[], monomers: string[][]): string {
  let helm = '';
  for (let i = 0; i < monomers.length; i++) {
    if (i > 0)
      helm += '|';

    helm += `PEPTIDE${i + 1}{`;

    for (let j = 0; j < monomers[i].length; j++) {
      if (j > 0)
        helm += '.';
      const symbol = monomers[i][j];
      helm += symbol.length > 1 ? `[${symbol}]` : symbol;
    }
    helm += `}`;
  }

  helm += '$';

  for (let i = 0; i < linkages.length; i++) {
    if (i > 0)
      helm += '|';
    helm += `PEPTIDE${linkages[i].fChain + 1},PEPTIDE${linkages[i].sChain + 1},`;

    helm += `${getInnerIdx(linkages[i].fMonomer - 1, monomers)[0] + 1}:R${linkages[i].fR}-`;
    helm += `${getInnerIdx(linkages[i].sMonomer - 1, monomers)[0] + 1}:R${linkages[i].sR}`;
  }

  helm += '$$$' + 'V2.0';
  return helm;
}

//homo and hetero dimers
export function handleDuplicated(sequence: string, rules: Rules): [Linkage[], string[]] {
  const mainFragments: string[] = [];
  const linkages: Linkage[] = [];
  const heterodimerCode = rules.heterodimerCode;
  const homodimerCode = rules.homodimerCode;

  //NOTICE: this works only with simple single heterodimers
  const heterodimeric = heterodimerCode !== null ? sequence.split(`(${rules.heterodimerCode!})`) : '';
  if (heterodimerCode !== null && heterodimeric.length > 1) {
    linkages.push({fChain: 0, sChain: 1, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
    mainFragments.push(heterodimeric[1].replaceAll('{', '').replaceAll('}', ''));
    mainFragments.push(heterodimeric[2].replaceAll('{', '').replaceAll('}', ''));
  } else {
    mainFragments.push(sequence);
  }

  //NOTICE: this works only with simple single dimers
  for (let i = 0; i < mainFragments.length; i++) {
    if (homodimerCode !== null && mainFragments[i].includes(`(${homodimerCode!})`)) {
      const idxSequence = mainFragments.length;

      linkages.push({fChain: i, sChain: idxSequence, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
      const rawDimer = mainFragments[i].replace(`(${homodimerCode!})`, '');
      const idx = rawDimer.indexOf('{');
      const linker = rawDimer.slice(0, idx);
      const body = rawDimer.replace(linker, '').replaceAll('{', '').replaceAll('}', '');

      mainFragments[i] = linker + body;
      mainFragments.push(body);
    }
  }

  for (let i = 0; i < mainFragments.length; i++) {
    if (homodimerCode !== null && mainFragments[i].includes(`(${homodimerCode!})`)) {
      const idxSequence = mainFragments.length;

      linkages.push({fChain: i, sChain: idxSequence, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
      const rawDimer = mainFragments[i].replace(`(${homodimerCode!})`, '');
      const idx = rawDimer.indexOf('{');
      const linker = rawDimer.slice(0, idx);
      const body = rawDimer.replace(linker, '').replaceAll('{', '').replaceAll('}', '');

      mainFragments[i] = linker + body;
      mainFragments.push(body);
    }
  }

  return [linkages, mainFragments];
}

export function handleLinkRules(mf: string[], monomers: string[][], linkages: Linkage[], rules: Rules): void {
  for (let i = 0; i < mf.length; i++) {
    const rawMonomers = mf[i].split('-');
    const linkedPositions = getLinkedPositions(rawMonomers, rules.linkRules);
    const [allPos1, allPos2, allAttaches1, allAttaches2] =
      getAllCycles(rules.linkRules, rawMonomers, linkedPositions);

    for (let j = 0; j < allPos1.length; j++) {
      linkages.push({
        fChain: i,
        sChain: i,
        fMonomer: allPos1[j],
        sMonomer: allPos2[j],
        fR: allAttaches1[j],
        sR: allAttaches2[j],
      });
    }

    monomers[i] = rawMonomers;
  }
}

export function handleReactionRules(monomers: string[][], linkages: Linkage[], rules: Rules): void {
  for (let i = 0; i < monomers.length; i++) {
    const linkedPositions = getLinkedPositions(monomers[i], rules.reactionRules);
    const [allPos1, allPos2, ruleN] = getAllReactants(rules.reactionRules, monomers[i], linkedPositions);

    if (allPos1.length >= 1) {
      linkages.push({
        fChain: i,
        sChain: monomers.length,
        fMonomer: allPos1[0],
        sMonomer: 1,
        fR: 3,
        sR: 1,
      });

      linkages.push({
        fChain: i,
        sChain: monomers.length,
        fMonomer: allPos2[0],
        sMonomer: 1,
        fR: 3,
        sR: 2,
      });

      monomers.push([rules.reactionRules[ruleN[0]].name]);
    }
  }
}

function reset(): [boolean, boolean, boolean, number, number] {
  return [false, false, false, -1, -1];
}

function find(monomers: string[], add: string, unlim: boolean, idx: number = 0,
  firstMonomers: string[] = [], secondMonomers: string[] = []) : [boolean, boolean, boolean, number, number] {
  let [firstFound, secondFound, firstIsFirst, firstIdx, secondIdx] = reset();
  for (let k = 0; k < monomers.length; k++) {
    if (monomers[k].includes(add)) {
      if (firstFound) {
        if (firstIsFirst && (unlim || monomers[k] == secondMonomers[idx] + add)) {
          secondFound = true;
          secondIdx = k;
          break;
        } else if (!firstIsFirst && (unlim || monomers[k] == firstMonomers[idx] + add)) {
          secondFound = true;
          secondIdx = k;
          break;
        } else {
          continue;
        }
      } else if (unlim) {
        firstFound = true;
        firstIsFirst = true;
        firstIdx = k;
      } else {
        if (monomers[k] == firstMonomers[idx] + add) {
          firstFound = true;
          firstIsFirst = true;
          firstIdx = k;
        } else if (monomers[k] == secondMonomers[idx] + add) {
          firstFound = true;
          firstIsFirst = unlim ? true : false;
          firstIdx = k;
        } else {
          continue;
        }
      }
    }
  }

  return [firstFound, secondFound, firstIsFirst, firstIdx, secondIdx];
}

function getLinkedPositions(monomers: string[], rules: RuleLink[] | RuleReaction []) : LinkedPosition[] {
  const result: LinkedPosition[] = [];

  for (let i = 0; i < rules.length; i++) {
    const add = `(${rules[i].code})`;

    const [firstMonomers, secondMonomers] = getMonomerPairs(rules[i]);

    if (firstMonomers.length > 0) {
      for (let j = 0; j < firstMonomers.length; j++) {
        const [firstFound, secondFound, firstIsFirst, firstIdx, secondIdx] =
          find(monomers, add, false, j, firstMonomers, secondMonomers);

        if (!(firstFound && secondFound))
          continue;
        else if (firstIsFirst)
          result.push({firstIdx, secondIdx, ruleIdx: i});
        else
          result.push({secondIdx, firstIdx, ruleIdx: i});
      }
    } else {
      const [firstFound, secondFound, firstIsFirst, firstIdx, secondIdx] = find(monomers, add, true);

      if (!(firstFound && secondFound))
        continue;
      else if (firstIsFirst)
        result.push({firstIdx, secondIdx, ruleIdx: i});
      else
        result.push({secondIdx, firstIdx, ruleIdx: i});
    }
  }

  return result;
}

function getAllCycles(rules: RuleLink[], monomers: string [], positions: LinkedPosition[]):
[number [], number [], number [], number []] {
  const allPos1: number [] = [];
  const allPos2: number [] = [];
  const allAttaches1: number [] = [];
  const allAttaches2: number [] = [];
  const count = positions.length;

  for (let i = 0; i < count; i++) {
    if (positions[i].firstIdx == -1)
      continue;

    const ruleNum = positions[i].ruleIdx;
    const code = rules[ruleNum].code;

    monomers[positions[i].firstIdx] = monomers[positions[i].firstIdx].replace(`(${code})`, '');
    monomers[positions[i].secondIdx] = monomers[positions[i].secondIdx].replace(`(${code})`, '');

    allPos1.push(positions[i].firstIdx + 1);
    allPos2.push(positions[i].secondIdx + 1);
    allAttaches1.push(rules[ruleNum].firstLinkingGroup);
    allAttaches2.push(rules[ruleNum].secondLinkingGroup);
  }

  return [allPos1, allPos2, allAttaches1, allAttaches2];
}

function getAllReactants(rules: RuleReaction[], monomers: string [], positions: LinkedPosition[]):
[number [], number [], number []] {
  const allPos1: number [] = [];
  const allPos2: number [] = [];
  const rule: number [] = [];
  const count = positions.length;

  for (let i = 0; i < count; i++) {
    if (positions[i].firstIdx == -1)
      continue;
    const fIdx = positions[i].firstIdx;
    const sIdx = positions[i].secondIdx;

    const ruleNum = positions[i].ruleIdx;
    const code = rules[ruleNum].code;
    monomers[fIdx] = monomers[fIdx].replace(`(${code})`, '') + `_${rules[ruleNum].name}`;
    monomers[sIdx] = monomers[sIdx].replace(`(${code})`, '') + `_${rules[ruleNum].name}`;

    allPos1.push(fIdx + 1);
    allPos2.push(sIdx + 1);
    rule.push(positions[i].ruleIdx);
  }

  return [allPos1, allPos2, rule];
}
