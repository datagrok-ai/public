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

class Chain {
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
  const fileSource = new DG.FileSource(RULES_PATH);
  const rulesRawFrames: DG.DataFrame[] = new Array<DG.DataFrame>(ruleFiles.length);

  for (let i = 0; i < ruleFiles.length; i++) {
    const rulesRaw = await fileSource.readAsText(ruleFiles[i].replace(RULES_PATH, ''));
    rulesRawFrames[i] = DG.DataFrame.fromCsv(rulesRaw);
  }

  const targetList = getHelms(molColumn.toList(), rulesRawFrames);
  const helmColName = df.columns.getUnusedName('transformed(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', helmColName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);

  const molCol = await grok.functions.call('Bio:getMolFromHelm', {
    'df': df,
    'helmCol': targetHelmCol,
    'chiralityEngine': chiralityEngine
  });


  molCol.name = df.columns.getUnusedName('molfile(' + molColumn.name + ')');
  molCol.semType = DG.SEMTYPE.MOLECULE;

  if (addHelm) {
    targetHelmCol.setTag('cell.renderer', 'helm');
    df.columns.add(targetHelmCol);
  }
  df.columns.add(molCol, true);
  await grok.data.detectSemanticTypes(df);
}
