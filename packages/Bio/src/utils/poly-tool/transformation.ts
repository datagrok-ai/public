import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {_package} from '../../package';
import {HELM_WRAPPER} from './const';
import {getMolColumnFromHelm} from '../helm-to-molfile';
import {ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

type ConnectionData = {
  allPos1: number[],
  allPos2: number[],
  allAttaches1: number[],
  allAttaches2: number[],
}

type Rule = {
  firstMonomer: string,
  secondMonomer: string,
  firstModification: string,
  secondModification: string,
  firstR: number,
  secondR: number
}

function addCommonTags(col: DG.Column):void {
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

  protected getLinkedPositions(helm: string, ruleCount: number): [number, number][] {
    const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
    const monomers = seq.split('.');
    const result:[number, number][] = new Array<[number, number]>(ruleCount);

    for (let i = 0; i < ruleCount; i++) {
      let firstFound = false;
      for (let j = 0; j < monomers.length; j++) {
        if (monomers[j].includes(`(${i + 1})`)) {
          if (firstFound) {
            result[i][1] = j;
            break;
          } else {
            firstFound = true;
            result[i] = [j, 0];
          }
        }
      }

      if (!firstFound)
        result[i] = [-1, -1];
    }


    return result;
  }

  protected getRules(rulesTable: DG.DataFrame): Rule[] {
    const ruleCount = rulesTable.rowCount;
    const rules: Rule[] = new Array<Rule>(ruleCount);
    const codeCol = rulesTable.columns.byName('code');
    const monomer1Col = rulesTable.columns.byName('monomer1');
    const monomer2Col = rulesTable.columns.byName('monomer2');
    const modification1Col = rulesTable.columns.byName('modification1');
    const modification2Col = rulesTable.columns.byName('modification2');
    const r1Col = rulesTable.columns.byName('R1');
    const r2Col = rulesTable.columns.byName('R2');

    for (let i = 0; i < ruleCount; i++) {
      rules[i] = {
        firstMonomer: monomer1Col.get(i),
        secondMonomer: monomer2Col.get(i),
        firstModification: modification1Col.get(i),
        secondModification: modification2Col.get(i),
        firstR: r1Col.get(i),
        secondR: r2Col.get(i),
      };
    }

    return rules;
  }

  protected getTransformedHelm(helm: string, rules: Rule[]): string {
    const ruleCount = rules.length;
    const positions = this.getLinkedPositions(helm, ruleCount);

    const allPos1: number [] = [];
    const allPos2: number [] = [];
    const allAttaches1: number [] = [];
    const allAttaches2: number [] = [];

    for (let i = 0; i < ruleCount; i++) {
      if (positions[i][0] == -1)
        continue;

      helm = helm.replaceAll(`(${i + 1})`, '');
      const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
      const monomers = seq.split('.');
      const firstMonomer = monomers[positions[i][0]].replace('[', '').replace(']', '');
      const secondMonomer = monomers[positions[i][1]].replace('[', '').replace(']', '');
      let attach1 = 0;
      let attach2 = 0;
      if (firstMonomer === rules[i].firstMonomer && secondMonomer === rules[i].secondMonomer) {
        monomers[positions[i][0]] = monomers[positions[i][0]].replace(firstMonomer, rules[i].firstModification);
        monomers[positions[i][1]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondModification);
        attach1 = rules[i].firstR;
        attach2 = rules[i].secondR;
      } else if (secondMonomer === rules[i].firstModification && firstMonomer === rules[i].secondModification) {
        monomers[positions[i][0]] = monomers[positions[i][1]].replace(secondMonomer, rules[i].secondModification);
        monomers[positions[i][1]] = rules[i].firstModification.replace(firstMonomer, rules[i].firstModification);
        attach1 = rules[i].secondR;
        attach2 = rules[i].firstR;
      } else {
        continue;
      }

      allPos1.push(positions[i][0] + 1);
      allPos2.push(positions[i][1] + 1);
      allAttaches1.push(attach1);
      allAttaches2.push(attach2);

      helm = HELM_WRAPPER.LEFT;
      for (let i = 0; i < monomers.length; i ++) {
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

  transform(rulesTable: DG.DataFrame): string[] {
    const rules = this.getRules(rulesTable);
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
  // return helm.replace(HELM_WRAPPER.RIGHT,
  //   `}$PEPTIDE1,PEPTIDE1,${
  //     source.monomerPosition
  //   }:R${
  //     source.attachmentPoint
  //   }-${
  //     target.monomerPosition
  //   }:R${
  //     target.attachmentPoint
  //   }${'$'.repeat(6)}`
  // );
}

export async function addTransformedColumn(
  molColumn: DG.Column<string>, rulesTable: DG.DataFrame, addHelm: boolean
): Promise<void> {
  const df = molColumn.dataFrame;
  const uh = UnitsHandler.getOrCreate(molColumn);
  const sourceHelmCol = uh.convert(NOTATION.HELM);
  const pt = PolymerTransformation.getInstance(sourceHelmCol);
  const targetList = pt.transform(rulesTable);
  const helmColName = df.columns.getUnusedName('transformed(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', helmColName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);

  const molCol = await getMolColumnFromHelm(df, targetHelmCol);
  molCol.name = df.columns.getUnusedName('molfile(' + molColumn.name + ')');

  if (addHelm) {
    targetHelmCol.setTag('cell.renderer', 'helm');
    df.columns.add(targetHelmCol);
  }
  df.columns.add(molCol, true);

  await grok.data.detectSemanticTypes(df);
}
