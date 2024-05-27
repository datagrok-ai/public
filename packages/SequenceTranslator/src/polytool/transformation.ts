import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

//import {getMolColumnFromHelm} from '../helm-to-molfile/utils';

export const RULES_PATH = 'System:AppData/Bio/polytool-rules/';
export const RULES_STORAGE_NAME = 'Polytool';

export const RULES_DIMER = '(#2)';
export const RULES_HETERODIMER = '($2)';

const enum HELM_WRAPPER {
  LEFT = 'PEPTIDE1{',
  RIGHT = '}$$$$',
}

type ConnectionData = {
  allPos1: number[],
  allPos2: number[],
  allAttaches1: number[],
  allAttaches2: number[],
}

type Rule = {
  code: number,
  firstMonomer: string,
  secondMonomer: string,
  firstModification: string,
  secondModification: string,
  firstR: number,
  secondR: number
}

function addCommonTags(col: DG.Column): void {
  col.setTag('quality', DG.SEMTYPE.MACROMOLECULE);
  col.setTag('aligned', ALIGNMENT.SEQ);
  col.setTag('alphabet', ALPHABET.PT);
}

class TransformationCommon {
  helmColumn: DG.Column;

  constructor(helmColumn: DG.Column<string>) {
    this.helmColumn = helmColumn;
  }

  protected getLinkedPositions(monomers: string[], rules: Rule[]): [number, number][] {
    const monomersNames = monomers.map((m) => { return m.replace('[', '').replace(']', ''); });
    const result: [number, number][] = new Array<[number, number]>(rules.length);

    for (let i = 0; i < rules.length; i++) {
      let firstFound = false;
      let secondFound = false;
      let firstIsFirst = false;
      let firstEntryIndex = -1;
      let secondEntryIndex = -1;
      const add = `(${rules[i].code})`;
      for (let j = 0; j < monomersNames.length; j++) {
        if (monomersNames[j].includes(add)) {
          if (firstFound) {
            if (firstIsFirst && monomersNames[j] == rules[i].secondMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else if (!firstIsFirst && monomersNames[j] == rules[i].firstMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else {
              continue;
            }
          } else {
            if (monomersNames[j] == rules[i].firstMonomer + add) {
              firstFound = true;
              firstIsFirst = true;
              firstEntryIndex = j;
            } else if (monomersNames[j] == rules[i].secondMonomer + add) {
              firstFound = true;
              firstIsFirst = false;
              firstEntryIndex = j;
            } else {
              continue;
            }
            //result[i] = [j, 0];
          }
        }
      }

      if (!(firstFound && secondFound))
        result[i] = [-1, -1];
      else if (firstIsFirst)
        result[i] = [firstEntryIndex, secondEntryIndex];
      else
        result[i] = [secondEntryIndex, firstEntryIndex];
    }


    return result;
  }

  protected getRules(rulesTables: DG.DataFrame[]): Rule[] {
    const ruleCount = rulesTables.map((df) => df.rowCount).reduce((a, b) => a + b);
    const rules: Rule[] = new Array<Rule>(ruleCount);

    let counter = 0;
    for (let i = 0; i < rulesTables.length; i++) {
      const codeCol = rulesTables[i].columns.byName('code');
      const monomer1Col = rulesTables[i].columns.byName('monomer1');
      const monomer2Col = rulesTables[i].columns.byName('monomer2');
      const modification1Col = rulesTables[i].columns.byName('modification1');
      const modification2Col = rulesTables[i].columns.byName('modification2');
      const r1Col = rulesTables[i].columns.byName('R1');
      const r2Col = rulesTables[i].columns.byName('R2');

      for (let j = 0; j < rulesTables[i].rowCount; j++, counter++) {
        rules[counter] = {
          code: codeCol.get(j),
          firstMonomer: monomer1Col.get(j),
          secondMonomer: monomer2Col.get(j),
          firstModification: modification1Col.get(j),
          secondModification: modification2Col.get(j),
          firstR: r1Col.get(j),
          secondR: r2Col.get(j),
        };
      }
    }

    return rules;
  }

  protected getDimeric(monomers: string[]): [number, number][] {
    const duplicates: [number, number][] = [];

    for (let i = 0; i < monomers.length; i++) {
      if (monomers[i].includes('(#2)')) {
        monomers[i] = monomers[i].replace('(#2)', '');
        const duplicateStart = i + 1;
        let duplicateFinish = 0;
        for (let j = duplicateStart + 1; j < monomers.length; j++) {
          if (monomers[j].includes('}')) {
            duplicateFinish = j; break;
          }
        }
        monomers[duplicateStart] = monomers[duplicateStart].replace('{', '');
        monomers[duplicateFinish] = monomers[duplicateFinish].replace('}', '');
        duplicates.push([duplicateStart, duplicateFinish]);
      }
    }

    return duplicates;
  }

  protected getHeteroDimeric(helm: string): [string[], [number, number, number][]] {
    const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
    const monomers = seq.split('.');
    const monomersHetero = new Array<string>(monomers.length);
    const heteros: [number, number, number][] = [];
    let delta = 0;

    for (let i = 0; i < monomers.length; i++) {
      if (monomers[i].includes('($2)')) {
        monomersHetero[i] = monomers[i].replace('($2)', '');
        const linkerIdx = i;

        while (!monomers[i].includes('}')) {
          i++;
          monomersHetero[i] = monomers[i].replace('{', '').replace('}', '');
        }
        const heteroFinish1 = i;
        monomers[i+1] = monomers[i+1].replace('($2)', '');
        while (!monomers[i].includes('}')) {
          i++;
          monomersHetero[i] = monomers[i].replace('{', '').replace('}', '');
        }
        const heteroFinish2 = i;
        heteros.push([linkerIdx, heteroFinish1, heteroFinish2]);
      } else {
        monomersHetero[i] = monomers[i];
      }
    }

    return [monomers, heteros];
  }

  protected getAllCycles(rules: Rule[], monomers: string [], positions: [number, number][]) :
  [string [], number [], number [], number [], number []] {
    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const allAttaches1: number [] = [];
    const allAttaches2: number [] = [];
    const ruleCount = rules.length;

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      const firstMonomer = monomers[positions[i][0]].replace('[', '').replace(']', '');
      const secondMonomer = monomers[positions[i][1]].replace('[', '').replace(']', '');

      monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstModification);
      monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondModification);

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      allAttaches1.push(rules[i].firstR);
      allAttaches2.push(rules[i].secondR);
    }

    return [monomers, allPos1, allPos2, allAttaches1, allAttaches2];
  }

  protected getHelmCycle(source: ConnectionData, polFirst: number, poSecond: number): string {
    let cycled = '';

    for (let i = 0; i < source.allPos1.length; i++) {
      if (i == 0)
        cycled += `PEPTIDE${polFirst},PEPTIDE${poSecond},`;
      else
        cycled += `|PEPTIDE${polFirst},PEPTIDE${poSecond},`;
      cycled += `${source.allPos1[i]}:R${source.allAttaches1[i]}-${source.allPos2[i]}:R${source.allAttaches2[i]}`;
    }
    return cycled;
  }

  protected getTransformedHelm(helm: string, rules: Rule[]): string {
    const [monomers, heteros] = this.getHeteroDimeric(helm);
    const duplicates = this.getDimeric(monomers);
    const positions = this.getLinkedPositions(monomers, rules);
    const [monomersCycled, allPos1, allPos2, allAttaches1, allAttaches2] =
      this.getAllCycles(rules, monomers, positions);

    helm = 'PEPTIDE1{';
    let heteroCounter = 0;
    for (let i = 0; i < monomersCycled.length; i++) {
      if (heteros[heteroCounter][1] + 1 === i) {
        i == heteros[heteroCounter][2];
        heteroCounter++;
      } else {
        helm += i != monomersCycled.length - 1 ? monomersCycled[i] + '.' : monomersCycled[i];
      }
    }

    helm += '}';

    const dimerCodes = new Array<string>(duplicates.length);
    const cycleCodes = new Array<string>(duplicates.length);
    for (let i = 0; i < duplicates.length; i++) {
      let helmAdd = `|PEPTIDE${i + 2}{`;

      const lengthAdd = duplicates[i][1] - duplicates[i][0];
      const allPosAdd1: number [] = [];
      const allPosAdd2: number [] = [];
      const allAttachesAdd1: number [] = [];
      const allAttachesAdd2: number [] = [];

      for (let j = 0; j <= lengthAdd; j ++) {
        const index = j + duplicates[i][0];
        helmAdd += j != lengthAdd ? monomersCycled[index] + '.' : monomersCycled[index];
      }

      helmAdd += '}';

      for (let j = 0; j < allPos1.length; j++) {
        if (allPos1[j] - 1 >= duplicates[i][0] && allPos1[j] - 1 <= duplicates[i][1]) {
          allPosAdd1.push(allPos1[j] - duplicates[i][0]);
          allPosAdd2.push(allPos2[j] - duplicates[i][0]);
          allAttachesAdd1.push(allAttaches1[j]);
          allAttachesAdd2.push(allAttaches1[j]);
        }
      }

      const addCyclysation = this.getHelmCycle({
        allPos1: allPosAdd1,
        allPos2: allPosAdd2,
        allAttaches1: allAttachesAdd1,
        allAttaches2: allAttachesAdd2}, i + 2, i + 2);

      dimerCodes[i] = helmAdd;
      cycleCodes[i] = this.getHelmCycle({
        allPos1: [duplicates[i][0]],
        allPos2: [1],
        allAttaches1: [1],
        allAttaches2: [1]}, 1, i + 2);

      if (addCyclysation !== '')
        cycleCodes.push(addCyclysation);
    }

    for (let i = 0; i < dimerCodes.length; i++)
      helm += dimerCodes[i];

    helm += '$';
    const mainCyclysation = this.getHelmCycle({allPos1, allPos2, allAttaches1, allAttaches2}, 1, 1);
    helm += mainCyclysation;
    for (let i = 0; i < cycleCodes.length; i++) {
      helm += (i === 0 && mainCyclysation === '') ? '' : '|';
      helm += cycleCodes[i];
    }
    helm += '$$$';
    return helm;
  }

  transform(rulesTables: DG.DataFrame[]): string[] {
    const rules = this.getRules(rulesTables);
    const resultList = this.helmColumn.toList().map((helm: string) => {
      return this.getTransformedHelm(helm, rules);
    });
    return resultList;
  }
}

// class PolymerTransformation {
//   private constructor() {}

//   static getInstance(molColumn: DG.Column<string>) {
//     return new TransformationCommon(molColumn);
//   }
// }

class Chain {
  //chains: Chain[];
  linkages: {fChain: number, sChain: number, fMonomer:number, sMonomer:number, fR:number, sR:number}[];
  monomers: string[][];

  constructor(sequence: string, rules: Rule[]) {
    const heterodimeric = sequence.split(RULES_HETERODIMER);
    const mainFragments: string[] = [];

    this.linkages = [];

    //NOTICE: this works only with simple single heterodimers
    if (heterodimeric.length > 1) {
      this.linkages.push({fChain: 0, sChain: 1, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
      mainFragments.push(heterodimeric[1].replaceAll('{', '').replaceAll('}', ''));
      mainFragments.push(heterodimeric[2].replaceAll('{', '').replaceAll('}', ''));
    } else {
      mainFragments.push(sequence);
    }

    //NOTICE: this works only with simple single dimers
    for (let i = 0; i < mainFragments.length; i++) {
      if (mainFragments[i].includes(RULES_DIMER)) {
        const idxSequence = mainFragments.length;

        this.linkages.push({fChain: i, sChain: idxSequence, fMonomer: 1, sMonomer: 1, fR: 1, sR: 1});
        const rawDimer = mainFragments[i].replace(RULES_DIMER, '');
        const idx = rawDimer.indexOf('{');
        const linker = rawDimer.slice(0, idx);
        const body = rawDimer.replace(linker, '').replaceAll('{', '').replaceAll('}', '');

        mainFragments[i] = linker + body;
        mainFragments.push(body);
      }
    }

    this.monomers = new Array<Array<string>>(mainFragments.length);

    //let monomerCounter = 0;
    for (let i = 0; i < mainFragments.length; i++) {
      const rawMonomers = mainFragments[i].split('-');
      const linkedPositions = this.getLinkedPositions(rawMonomers, rules);
      const [monomersCycled, allPos1, allPos2, allAttaches1, allAttaches2] =
      this.getAllCycles(rules, rawMonomers, linkedPositions);

      const monomersReady = new Array<string>(monomersCycled.length);
      for (let j = 0; j < monomersCycled.length; j++)
        monomersReady[j] = `[${monomersCycled[j]}]`;

      for (let j = 0; j < allPos1.length; j++) {
        this.linkages.push({
          fChain: i,
          sChain: i,
          fMonomer: allPos1[j],
          sMonomer: allPos2[j],
          fR: allAttaches1[j],
          sR: allAttaches2[j],
        });
      }

      this.monomers[i] = monomersReady;
      //monomerCounter += monomersCycled.length;
    }
  }

  getHelm(): string {
    let helm = '';
    for (let i = 0; i < this.monomers.length; i++) {
      if (i > 0)
        helm += '|';

      helm += `PEPTIDE${i + 1}{`;

      for (let j = 0; j < this.monomers[i].length; j++) {
        if (j > 0)
          helm += '.';
        helm += this.monomers[i][j];
      }
      helm += `}`;
    }

    helm += '$';

    for (let i = 0; i < this.linkages.length; i++) {
      if (i > 0)
        helm += '|';
      helm += `PEPTIDE${this.linkages[i].fChain + 1},PEPTIDE${this.linkages[i].sChain + 1},`;
      helm += `${this.linkages[i].fMonomer}:R${this.linkages[i].fR}-`;
      helm += `${this.linkages[i].sMonomer}:R${this.linkages[i].sR}`;
    }

    helm += '$$$';
    return helm;
  }

  protected getLinkedPositions(monomers: string[], rules: Rule[]): [number, number][] {
    const result: [number, number][] = new Array<[number, number]>(rules.length);

    for (let i = 0; i < rules.length; i++) {
      let firstFound = false;
      let secondFound = false;
      let firstIsFirst = false;
      let firstEntryIndex = -1;
      let secondEntryIndex = -1;
      const add = `(${rules[i].code})`;
      for (let j = 0; j < monomers.length; j++) {
        if (monomers[j].includes(add)) {
          if (firstFound) {
            if (firstIsFirst && monomers[j] == rules[i].secondMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else if (!firstIsFirst && monomers[j] == rules[i].firstMonomer + add) {
              secondFound = true;
              secondEntryIndex = j;
              break;
            } else {
              continue;
            }
          } else {
            if (monomers[j] == rules[i].firstMonomer + add) {
              firstFound = true;
              firstIsFirst = true;
              firstEntryIndex = j;
            } else if (monomers[j] == rules[i].secondMonomer + add) {
              firstFound = true;
              firstIsFirst = false;
              firstEntryIndex = j;
            } else {
              continue;
            }
          }
        }
      }

      if (!(firstFound && secondFound))
        result[i] = [-1, -1];
      else if (firstIsFirst)
        result[i] = [firstEntryIndex, secondEntryIndex];
      else
        result[i] = [secondEntryIndex, firstEntryIndex];
    }


    return result;
  }

  protected getAllCycles(rules: Rule[], monomers: string [], positions: [number, number][]) :
  [string [], number [], number [], number [], number []] {
    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const allAttaches1: number [] = [];
    const allAttaches2: number [] = [];
    const ruleCount = rules.length;

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      const firstMonomer = monomers[positions[i][0]];
      const secondMonomer = monomers[positions[i][1]];

      monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstModification);
      monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondModification);

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      allAttaches1.push(rules[i].firstR);
      allAttaches2.push(rules[i].secondR);
    }

    return [monomers, allPos1, allPos2, allAttaches1, allAttaches2];
  }
}

function getRules(rulesTables: DG.DataFrame[]): Rule[] {
  const ruleCount = rulesTables.map((df) => df.rowCount).reduce((a, b) => a + b);
  const rules: Rule[] = new Array<Rule>(ruleCount);

  let counter = 0;
  for (let i = 0; i < rulesTables.length; i++) {
    const codeCol = rulesTables[i].columns.byName('code');
    const monomer1Col = rulesTables[i].columns.byName('monomer1');
    const monomer2Col = rulesTables[i].columns.byName('monomer2');
    const modification1Col = rulesTables[i].columns.byName('modification1');
    const modification2Col = rulesTables[i].columns.byName('modification2');
    const r1Col = rulesTables[i].columns.byName('R1');
    const r2Col = rulesTables[i].columns.byName('R2');

    for (let j = 0; j < rulesTables[i].rowCount; j++, counter++) {
      rules[counter] = {
        code: codeCol.get(j),
        firstMonomer: monomer1Col.get(j),
        secondMonomer: monomer2Col.get(j),
        firstModification: modification1Col.get(j),
        secondModification: modification2Col.get(j),
        firstR: r1Col.get(j),
        secondR: r2Col.get(j),
      };
    }
  }

  return rules;
}


function getHelms(sequences: string[], rulesTables: DG.DataFrame[]): string[] {
  const rules = getRules(rulesTables);
  const helms = new Array<string>(sequences.length);
  for (let i = 0; i < sequences.length; i++) {
    const chain = new Chain(sequences[i], rules);
    helms[i] = chain.getHelm();
  }


  return helms;
}

export async function addTransformedColumn(
  molColumn: DG.Column<string>, addHelm: boolean, ruleFiles: string[], chiralityEngine?: boolean
): Promise<void> {
  const df = molColumn.dataFrame;
  //const sh = SeqHandler.forColumn(molColumn);
  //const sourceHelmCol = sh.convert(NOTATION.HELM);
  //const pt = PolymerTransformation.getInstance(sourceHelmCol);
  const fileSource = new DG.FileSource(RULES_PATH);
  const rulesRawFrames: DG.DataFrame[] = new Array<DG.DataFrame>(ruleFiles.length);

  for (let i = 0; i < ruleFiles.length; i++) {
    const rulesRaw = await fileSource.readAsText(ruleFiles[i].replace(RULES_PATH, ''));
    rulesRawFrames[i] = DG.DataFrame.fromCsv(rulesRaw);
  }

  const targetList = getHelms(molColumn.toList(), rulesRawFrames);
  //const targetList = pt.transform(rulesRawFrames);
  const helmColName = df.columns.getUnusedName('transformed(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', helmColName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);

  //const molCol = await getMolColumnFromHelm(df, targetHelmCol, chiralityEngine);
  const molCol = await grok.functions.call('Bio:getMolFromHelm', {
    'df': df,
    'helmCol': targetHelmCol,
    'chiralityEngine': chiralityEngine
  });


  molCol.name = df.columns.getUnusedName('molfile(' + molColumn.name + ')');
  molCol.semType = DG.SEMTYPE.MOLECULE;

  if (addHelm) {
    targetHelmCol.setTag('cell.renderer', 'helm');
    //targetHelmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    df.columns.add(targetHelmCol);
  }
  df.columns.add(molCol, true);
  await grok.data.detectSemanticTypes(df);
}
