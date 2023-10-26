import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {_package} from '../../package';
import {addCommonTags} from './utils';
import {HELM_WRAPPER, ALL_MONOMERS, CYCLIZATION_TYPE} from './const';
import {MetaData, ConnectionData} from './types';
import {getMolColumnFromHelm} from '../helm-to-molfile';

abstract class TransformationBase {
  constructor(helmColumn: DG.Column<string>, meta: MetaData) {
    this.helmColumn = helmColumn;
    this.leftTerminal = meta.leftTerminal;
    this.rightTerminal = meta.rightTerminal;
  }

  protected helmColumn: DG.Column<string>;
  protected leftTerminal?: string;
  protected rightTerminal?: string;

  protected abstract hasTerminals(helm: string): boolean;
  protected abstract getTransformedHelm(helm: string): string;
  protected abstract getLinkedPositions(helm: string): [number, number];

  transform(): string[] {
    const resultList = this.helmColumn.toList().map((helm: string) => {
      if (this.hasTerminals(helm))
        return this.getTransformedHelm(helm);
      return helm;
    });
    return resultList;
  }
}

class TransformationNCys extends TransformationBase {
  constructor(helmColumn: DG.Column<string>, meta: MetaData) {
    super(helmColumn, meta);
  }

  protected hasTerminals(helm: string): boolean {
    // if (! helm.includes(this.rightTerminal + HELM_WRAPPER.RIGHT))
    //   return false;
    // if (this.leftTerminal === ALL_MONOMERS)
    //   return true;
    // return helm.includes(HELM_WRAPPER.LEFT + this.leftTerminal);
    const positions = this.getLinkedPositions(helm);
    return positions.every((el) => el > 0);
  }

  protected getLinkedPositions(helm: string): [number, number] {
    const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
    const monomers = seq.split('.');
    const start = 0;
    const end = monomers.findIndex((el, idx) => el === this.rightTerminal && idx > start);
    return [start + 1, end + 1];
  }

  protected getTransformedHelm(helm: string): string {
    const positions = this.getLinkedPositions(helm);
    const source = {monomerPosition: positions[0], attachmentPoint: 1};
    const target = {monomerPosition: positions[1], attachmentPoint: 3};
    return getHelmCycle(helm, source, target);
  }
}

class TransformationNO extends TransformationBase {
  constructor(helmColumn: DG.Column<string>, meta: MetaData) {
    super(helmColumn, meta);
  }

  protected hasTerminals(helm: string): boolean {
    if (this.leftTerminal === ALL_MONOMERS || this.rightTerminal === ALL_MONOMERS)
      return true;
    return helm.includes(HELM_WRAPPER.LEFT + this.leftTerminal) &&
      helm.includes(this.rightTerminal + HELM_WRAPPER.RIGHT);
  }

  protected getLinkedPositions(helm: string): [number, number] {
    return [1, getNumberOfMonomers(helm)];
  }

  protected getTransformedHelm(helm: string): string {
    const positions = this.getLinkedPositions(helm);
    const source = {monomerPosition: positions[0], attachmentPoint: 1};
    const target = {monomerPosition: positions[1], attachmentPoint: 2};
    return getHelmCycle(helm, source, target);
  }
}

class TransformationR3 extends TransformationBase {
  constructor(helmColumn: DG.Column<string>, meta: MetaData) {
    super(helmColumn, meta);
  }

  protected hasTerminals(helm: string): boolean {
    if (this.leftTerminal === ALL_MONOMERS || this.rightTerminal === ALL_MONOMERS)
      return true;
    const positions = this.getLinkedPositions(helm);
    return positions.every((el) => el > 0);
  }

  protected getLinkedPositions(helm: string): [number, number] {
    const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
    const monomers = seq.split('.');
    const start = monomers.findIndex((el) => el === this.leftTerminal);
    const end = monomers.findIndex((el, idx) => el === this.rightTerminal && idx > start);
    return [start + 1, end + 1];
  }

  protected getTransformedHelm(helm: string): string {
    const positions = this.getLinkedPositions(helm);
    const source = {monomerPosition: positions[0], attachmentPoint: 3};
    const target = {monomerPosition: positions[1], attachmentPoint: 3};
    return getHelmCycle(helm, source, target);
  }
}

class PolymerTransformation {
  private constructor() {}

  static getInstance(molColumn: DG.Column<string>, meta: MetaData): TransformationBase {
    const cyclizationType = meta.cyclizationType;
    return (cyclizationType === CYCLIZATION_TYPE.R3) ? new TransformationR3(molColumn, meta) :
      (cyclizationType === CYCLIZATION_TYPE.NCys) ? new TransformationNCys(molColumn, meta) :
        new TransformationNO(molColumn, meta);
  }
}

function getNumberOfMonomers(helm: string): number {
  const seq = helm.replace(HELM_WRAPPER.LEFT, '').replace(HELM_WRAPPER.RIGHT, '');
  return seq.split('.').length;
}

function getHelmCycle(helm: string, source: ConnectionData, target: ConnectionData): string {
  return helm.replace(HELM_WRAPPER.RIGHT,
    `}$PEPTIDE1,PEPTIDE1,${
      source.monomerPosition
    }:R${
      source.attachmentPoint
    }-${
      target.monomerPosition
    }:R${
      target.attachmentPoint
    }${'$'.repeat(6)}`
  );
}

export async function addTransformedColumn(
  molColumn: DG.Column<string>, meta: MetaData, addHelm: boolean
): Promise<void> {
  const df = molColumn.dataFrame;
  const uh = UnitsHandler.getOrCreate(molColumn);
  const sourceHelmCol = uh.convert(NOTATION.HELM);
  const pt = PolymerTransformation.getInstance(sourceHelmCol, meta);
  const targetList = pt.transform();
  const helmColName = df.columns.getUnusedName(`${meta.transformationType}(` + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', helmColName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);

  const molCol = await getMolColumnFromHelm(df, targetHelmCol);
  molCol.name = df.columns.getUnusedName(`${meta.transformationType}_molfile(` + molColumn.name + ')');

  if (addHelm) {
    targetHelmCol.setTag('cell.renderer', 'helm');
    df.columns.add(targetHelmCol);
  }
  df.columns.add(molCol, true);

  await grok.data.detectSemanticTypes(df);
}
