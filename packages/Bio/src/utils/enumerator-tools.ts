
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibHelper} from '../utils/monomer-lib';
import {_package} from '../package';
import * as rxjs from 'rxjs';

const LEFT_HELM_WRAPPER = 'PEPTIDE1{';
const RIGHT_HELM_WRAPPER = '}$$$$';
const ALL_MONOMERS = '<All>';

type MetaData = {
  leftTerminal: string,
  rightTerminal: string,
  transformationType: TRANSFORMATION_TYPE,
  cyclizationType: CYCLIZATION_TYPE,
}

type ConnectionData = {
  monomerPosition: number,
  attachmentPoint: number,
}

const enum TRANSFORMATION_TYPE {
  CYCLIZATION = 'Cyclization',
}

const enum CYCLIZATION_TYPE {
  NO = 'N-O',
  NCys = 'N-Cys',
  R3 = 'R3-R3',
}

function addCommonTags(col: DG.Column):void {
  col.setTag('quality', DG.SEMTYPE.MACROMOLECULE);
  col.setTag('aligned', ALIGNMENT.SEQ);
  col.setTag('alphabet', ALPHABET.PT);
}

export function _setPeptideColumn(col: DG.Column): void {
  addCommonTags(col);
  col.setTag('units', NOTATION.SEPARATOR);
  col.setTag('separator', '-');
  // col.setTag('cell.renderer', 'sequence');
}

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
    if (! helm.includes(this.rightTerminal + RIGHT_HELM_WRAPPER))
      return false;
    if (this.leftTerminal === ALL_MONOMERS)
      return true;
    return helm.includes(LEFT_HELM_WRAPPER + this.leftTerminal);
  }

  protected getLinkedPositions(helm: string): [number, number] {
    return [1, getNumberOfMonomers(helm)];
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
    return helm.includes(LEFT_HELM_WRAPPER + this.leftTerminal) &&
      helm.includes(this.rightTerminal + RIGHT_HELM_WRAPPER);
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
    const seq = helm.replace(LEFT_HELM_WRAPPER, '').replace(RIGHT_HELM_WRAPPER, '');
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
  const seq = helm.replace(LEFT_HELM_WRAPPER, '').replace(RIGHT_HELM_WRAPPER, '');
  return seq.split('.').length;
}

function getHelmCycle(helm: string, source: ConnectionData, target: ConnectionData): string {
  return helm.replace(RIGHT_HELM_WRAPPER,
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

async function addTransformedColumn(
  molColumn: DG.Column<string>, meta: MetaData
): Promise<void> {
  const df = molColumn.dataFrame;
  const uh = UnitsHandler.getOrCreate(molColumn);
  const sourceHelmCol = uh.convert(NOTATION.HELM);
  const pt = PolymerTransformation.getInstance(sourceHelmCol, meta);
  const targetList = pt.transform();
  const colName = df.columns.getUnusedName(`${meta.transformationType}(` + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', colName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);
  targetHelmCol.setTag('cell.renderer', 'helm');

  df.columns.add(targetHelmCol);
  await grok.data.detectSemanticTypes(df);
}

export function _getEnumeratorWidget(molColumn: DG.Column): DG.Widget {
  function getMonomerList(cyclizationType: CYCLIZATION_TYPE): string[] {
    if (cyclizationType === cyclizationTypes[0]) {
      return [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE)
      );
    }
    if (cyclizationType === cyclizationTypes[1]) {
      return [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByRGroup(3, HELM_POLYMER_TYPE.PEPTIDE)
      );
    }
    return ['C'];
  }

  function updateMonomerList(): void {
    if (cyclizationTypeChoice.value === CYCLIZATION_TYPE.NCys) {
      monomerList1 = getMonomerList(CYCLIZATION_TYPE.NO);
      monomerList2 = getMonomerList(CYCLIZATION_TYPE.NCys);
    } else {
      monomerList1 = getMonomerList(cyclizationTypeChoice.value as CYCLIZATION_TYPE);
      monomerList2 = [...monomerList1];
    }

    leftTerminalChoice = ui.choiceInput(
      'R1:', monomerList1[0], monomerList1, () => { onRGroupValueChange.next(); }
    );
    rightTerminalChoice = ui.choiceInput('R2:', monomerList2[0], monomerList2, () => { onRGroupValueChange.next(); });
    onRGroupValueChange.next();
    ui.empty(terminalControls);
    [leftTerminalChoice, rightTerminalChoice].forEach((el) => { terminalControls.appendChild(el.root); });
  }

  const onCyclizationChoice = new rxjs.Subject<string>();
  const onRGroupValueChange = new rxjs.Subject<string>();
  onCyclizationChoice.subscribe(() => {
    meta.cyclizationType = cyclizationTypeChoice.value!;
    updateMonomerList();
  });
  onRGroupValueChange.subscribe(() => {
    meta.rightTerminal = rightTerminalChoice.value!;
    meta.leftTerminal = leftTerminalChoice.value!;
  });


  const meta = {} as MetaData;
  const transformations = [TRANSFORMATION_TYPE.CYCLIZATION];
  const transformationChoice = ui.choiceInput(
    'Modification', transformations[0], transformations, () => meta.transformationType = transformationChoice.value!
  );

  const cyclizationTypes = [CYCLIZATION_TYPE.NO, CYCLIZATION_TYPE.R3, CYCLIZATION_TYPE.NCys];
  const cyclizationTypeChoice = ui.choiceInput(
    'Type', cyclizationTypes[0], cyclizationTypes, () => { onCyclizationChoice.next(); }
  );

  const monomerLib = MonomerLibHelper.instance.getBioLib();
  let monomerList1: string[] = [];
  let monomerList2: string[] = [];
  let leftTerminalChoice = ui.choiceInput(
    'R1:', monomerList1[0], monomerList1, () => {
      meta.leftTerminal = leftTerminalChoice.value!;
    }
  );
  let rightTerminalChoice = ui.choiceInput('R2:', monomerList2[0], monomerList2, () => {
    meta.rightTerminal = rightTerminalChoice.value!;
  });
  const terminalControls = ui.divV([leftTerminalChoice.root, rightTerminalChoice.root]);

  function updateMeta() {
    meta.cyclizationType = cyclizationTypeChoice.value!;
    meta.leftTerminal = leftTerminalChoice.value!;
    meta.rightTerminal = rightTerminalChoice.value!;
    meta.transformationType = transformationChoice.value!;
  }

  updateMonomerList();

  updateMeta();

  const btn = ui.bigButton('Run', async () =>
    addTransformedColumn(molColumn, meta)
  );

  const div = ui.div([
    transformationChoice,
    cyclizationTypeChoice,
    terminalControls,
    btn
  ]);

  return new DG.Widget(div);
}
