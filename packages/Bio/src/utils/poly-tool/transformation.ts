import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {HELM_WRAPPER} from './const';
import {getMolColumnFromHelm} from '../helm-to-molfile';

export const RULES_PATH = 'System:AppData/Bio/polytool-rules/';
export const RULES_STORAGE_NAME = 'Polytool';


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

  protected hasTerminals(helm: string): boolean {
    let isLinkable = false;
    if (helm.includes('(1)'))
      isLinkable = true;
    if (helm.includes('(2)'))
      isLinkable = true;

    return isLinkable;
  }

  protected getLinkedPositions(helm: string, rules: Rule[]): [number, number][] {
    const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
    const monomers = seq.split('.').map((m) => { return m.replace('[', '').replace(']', ''); });
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
            //result[i][1] = j;
            // secondFound = true;
            // break;
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

  protected getTransformedHelm(helm: string, rules: Rule[]): string {
    const ruleCount = rules.length;
    const positions = this.getLinkedPositions(helm, rules);

    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const allAttaches1: number [] = [];
    const allAttaches2: number [] = [];

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      //helm = helm.replaceAll(`(${i + 1})`, '');
      const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
      const monomers = seq.split('.');
      const firstMonomer = monomers[positions[i][0]].replace('[', '').replace(']', '');
      const secondMonomer = monomers[positions[i][1]].replace('[', '').replace(']', '');

      monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstModification);
      monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondModification);

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      allAttaches1.push(rules[i].firstR);
      allAttaches2.push(rules[i].secondR);

      helm = HELM_WRAPPER.LEFT;
      for (let i = 0; i < monomers.length; i++) {
        if (i != monomers.length - 1)
          helm = helm + monomers[i] + '.';
        else
          helm = helm + monomers[i];
      }
      helm = helm + HELM_WRAPPER.RIGHT;
    }


    const cycledHelm = getHelmCycle(helm, {allPos1, allPos2, allAttaches1, allAttaches2});

    return cycledHelm;
  }

  transform(rulesTables: DG.DataFrame[]): string[] {
    const rules = this.getRules(rulesTables);
    const resultList = this.helmColumn.toList().map((helm: string) => {
      if (this.hasTerminals(helm))
        return this.getTransformedHelm(helm, rules);

      console.log(helm);
      return helm;
    });
    return resultList;
  }
}

class PolymerTransformation {
  private constructor() {}

  static getInstance(molColumn: DG.Column<string>) {
    return new TransformationCommon(molColumn);
  }
}

function getHelmCycle(helm: string, source: ConnectionData): string {
  let cycled = helm.replace(HELM_WRAPPER.RIGHT, '}$');

  for (let i = 0; i < source.allPos1.length; i++) {
    if (i == 0)
      cycled += 'PEPTIDE1,PEPTIDE1,';
    else
      cycled += '|PEPTIDE1,PEPTIDE1,';
    cycled += `${source.allPos1[i]}:R${source.allAttaches1[i]}-${source.allPos2[i]}:R${source.allAttaches2[i]}`;
  }

  cycled += '$$$';
  return cycled;
}

export async function addTransformedColumn(
  molColumn: DG.Column<string>, addHelm: boolean, ruleFiles: string[], chiralityEngine?: boolean
): Promise<void> {
  const df = molColumn.dataFrame;
  const sh = SeqHandler.forColumn(molColumn);
  const sourceHelmCol = sh.convert(NOTATION.HELM);
  const pt = PolymerTransformation.getInstance(sourceHelmCol);
  const fileSource = new DG.FileSource(RULES_PATH);

  const rulesRawFrames: DG.DataFrame[] = new Array<DG.DataFrame>(ruleFiles.length);

  for (let i = 0; i < ruleFiles.length; i++) {
    const rulesRaw = await fileSource.readAsText(ruleFiles[i].replace(RULES_PATH, ''));
    rulesRawFrames[i] = DG.DataFrame.fromCsv(rulesRaw);
  }

  const targetList = pt.transform(rulesRawFrames);
  const helmColName = df.columns.getUnusedName('transformed(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', helmColName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);

  const molCol = await getMolColumnFromHelm(df, targetHelmCol, chiralityEngine);
  molCol.name = df.columns.getUnusedName('molfile(' + molColumn.name + ')');

  if (addHelm) {
    targetHelmCol.setTag('cell.renderer', 'helm');
    df.columns.add(targetHelmCol);
  }
  df.columns.add(molCol, true);

  await grok.data.detectSemanticTypes(df);
}
