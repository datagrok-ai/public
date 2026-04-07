import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { IPdbHelper, getPdbHelper } from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import { MolfileHandler } from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import '../../css/bionemo.css';

// Constants
export const CONSTANTS = {
  POSES_COLUMN_NAME: 'poses',
  CONFIDENCE_COLUMN_NAME: 'confidence',
  VIRTUAL_POSES_COLUMN_NAME: '~posesJson',
  TARGET_PATH: 'System:AppData/Bionemo/targets',
  GRID_ROW_HEIGHT: 100,
  POSES_COLUMN_WIDTH: 100,
  CONFIDENCE_COLUMN_WIDTH: 150
};

export interface PosesJson {
  ligand_positions: string[];
  position_confidence: number[];
  receptor: string;
}

export class DiffDockModel {
  private df: DG.DataFrame;
  private ligands: DG.Column;
  private target: string;
  private targetName: string;
  private poses: number;
  private currentViewer: DG.DockNode | null = null;
  public posesColumnName!: string;
  public virtualPosesColumnName!: string;

  constructor(df: DG.DataFrame, ligands: DG.Column, target: string, targetName: string, poses: number) {
    this.df = df;
    this.ligands = ligands;
    this.target = target;
    this.targetName = targetName;
    this.poses = poses;
  }

  async run(): Promise<DG.DataFrame> {
    const result = await this.calculatePoses();
    this.subscribeToCurrentCellChanged();
    return result;
  }

  private async calculatePoses(): Promise<DG.DataFrame> {
    const posesColumn = this.createColumn(DG.TYPE.STRING, CONSTANTS.POSES_COLUMN_NAME, this.ligands.length);
    const confidenceColumn = this.createColumn(DG.TYPE.FLOAT, CONSTANTS.CONFIDENCE_COLUMN_NAME, this.ligands.length);
    const virtualPosesColumn = this.createColumn(DG.TYPE.STRING, `${CONSTANTS.VIRTUAL_POSES_COLUMN_NAME}_${this.targetName}_${this.poses}`, this.ligands.length);
    const grid = grok.shell.getTableView(this.df.name).grid;

    for (let i = 0; i < posesColumn.length; ++i) {
      const posesJson = await this.getPosesJson(this.ligands.get(i), this.ligands.getTag(DG.TAGS.UNITS));
      const pdbHelper: IPdbHelper = await getPdbHelper();
      const { bestId, bestPose, confidence } = this.findBestPose(posesJson);

      posesColumn.set(i, await pdbHelper.molToPdb(bestPose));
      posesColumn.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
      posesColumn.setTag('docking.role', 'ligand');

      confidenceColumn.set(i, confidence);
      confidenceColumn.setTag(DG.TAGS.FORMAT, '0.00');
      confidenceColumn.meta.colors.setLinear([DG.Color.green, DG.Color.red]);
      grid.sort([confidenceColumn]);

      posesJson.receptor = this.target;
      virtualPosesColumn.set(i, JSON.stringify(posesJson));
    }

    this.posesColumnName = posesColumn.name;
    this.virtualPosesColumnName = virtualPosesColumn.name;

    grid.onCellRender.subscribe(() => this.configureGrid(grid));
    
    const resultDf = DG.DataFrame.fromColumns([posesColumn, confidenceColumn, virtualPosesColumn]);
    await grok.data.detectSemanticTypes(resultDf);
    return resultDf;
  }

  private configureGrid(grid: DG.Grid) {
    grid.setOptions({ 'rowHeight': CONSTANTS.GRID_ROW_HEIGHT });
    grid.col(CONSTANTS.POSES_COLUMN_NAME)!.width = CONSTANTS.POSES_COLUMN_WIDTH;
    grid.col(CONSTANTS.CONFIDENCE_COLUMN_NAME)!.width = CONSTANTS.CONFIDENCE_COLUMN_WIDTH;
    grid.col(CONSTANTS.CONFIDENCE_COLUMN_NAME)!.isTextColorCoded = true;
    grid.invalidate();
  }

  public async getPosesJson(ligand: string, units: string): Promise<PosesJson> {
    const isSmiles = units === DG.UNITS.Molecule.SMILES;
    const ligandValue = isSmiles ? DG.chem.convert(ligand, DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock) : ligand;
    const sdf = MolfileHandler.getInstance(ligandValue).z.every((coord) => coord === 0) ?
      (await grok.functions.call('Chem:SmilesTo3DCoordinates', { molecule: ligandValue })).replaceAll('\\n', '\n') : ligandValue;
    const jsonText = await grok.functions.call('Bionemo:diffDockModelScript', { ligand: sdf, target: this.target, poses: this.poses });
    return JSON.parse(jsonText);
  }

  public createColumn(type: DG.ColumnType, name: string, length: number): DG.Column {
    const unusedName = this.df.columns.getUnusedName(name);
    return DG.Column.fromType(type, unusedName, length);
  }

  public findBestPose(poses: PosesJson) {
    const positions = poses.position_confidence;
    const minValue = Math.min(...positions);
    const idx = positions.findIndex((value) => value === minValue);
    return {
      bestId: idx,
      bestPose: poses.ligand_positions[idx],
      confidence: minValue,
    };
  }

  private createPoseItems(posesJson: PosesJson): string[] {
    return posesJson.ligand_positions.map((_, index) =>
      `Pose ${index + 1} - Score: ${posesJson.position_confidence[index].toFixed(3)}`
    );
  }

  public async createCombinedControl(poses: PosesJson, bestPose: string, bestId: number, fullHeight: boolean = true): Promise<HTMLDivElement> {
    const molstarViewer = await this.createMolstarViewer(bestPose);
    molstarViewer.root.classList.add('bsv-container-info-panel');
    molstarViewer.root.style.alignSelf = 'center';
    if (fullHeight)
      molstarViewer.root.style.height = '100%';

    const items = this.createPoseItems(poses);
    const molstarProps = molstarViewer.getProperties();
    const ligandValueProp = molstarProps.find((p: DG.Property) => p.name === 'ligandValue');
    const ligandColProp = molstarProps.find((p: DG.Property) => p.name === 'ligandColumnName');

    const comboPopup = ui.comboPopup(items[bestId], items, (item: string) => {
      comboPopup.getElementsByTagName('span')[0].textContent = item;
      const idx = parseInt(item.split(' ')[1]) - 1;
      const ligandObject = {
        value: poses.ligand_positions[idx],
        semType: DG.SEMTYPE.MOLECULE,
        units: DG.UNITS.Molecule.MOLBLOCK,
      };
      molstarViewer.apply({ 'ligandValue': JSON.stringify(ligandObject) });
      molstarViewer.onPropertyChanged(ligandValueProp!);
      molstarViewer.onPropertyChanged(ligandColProp!);
    });

    comboPopup.classList.add('bionemo-combo-popup');
    const container = ui.divV([comboPopup, molstarViewer.root]);
    container.style.height = '100%';
    return container;
  }

  private async createMolstarViewer(bestPose: string): Promise<DG.Widget> {
    const ligandObject = {
      value: bestPose,
      semType: DG.SEMTYPE.MOLECULE,
      units: DG.UNITS.Molecule.MOLBLOCK,
    };

    return await this.ligands.dataFrame.plot.fromType('Biostructure', {
      pdb: this.target,
      ligandColumnName: this.ligands.name,
      ligandValue: JSON.stringify(ligandObject),
      zoom: true,
    });
  }

  private async handleCurrentCellChanged() {
    const currentCell = this.df.currentCell;
    const posesColumn = this.df.columns.byName(this.posesColumnName);
    if (currentCell && currentCell.column === posesColumn) {
      const currentRow = this.df.currentRowIdx;
      const virtualPosesColumn = this.df.columns.byName(this.virtualPosesColumnName);
      const poses = (JSON.parse(virtualPosesColumn.get(currentRow)) as PosesJson);
      const { bestId, bestPose, confidence } = this.findBestPose(poses);
      const combinedControl = await this.createCombinedControl(poses, bestPose, bestId);

      const view = grok.shell.getTableView(this.df.name);
      if (this.currentViewer)
        view.dockManager.close(this.currentViewer);

      this.currentViewer = view.dockManager.dock(ui.divV([combinedControl]), DG.DOCK_TYPE.RIGHT, null, 'Mol*', 0.35);
    }
  }

  public subscribeToCurrentCellChanged() {
    this.df.onCurrentCellChanged.subscribe(() => this.handleCurrentCellChanged());
  }
}