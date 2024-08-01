import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {IPdbHelper, getPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import '../css/bionemo.css';

interface PosesJson {
  ligand_positions: string[];
  position_confidence: number[];
  receptor: string;
}
  
export const CONSTANTS = {
  POSES_COLUMN_NAME: 'poses',
  CONFIDENCE_COLUMN_NAME: 'confidence',
  VIRTUAL_POSES_COLUMN_NAME: '~popesJson'
};
  
export class DiffDockModel {
  private df: DG.DataFrame;
  private ligands: DG.Column;
  private target: string;
  private poses: number;
  public posesColumn!: DG.Column;
  public virtualPosesColumn!: DG.Column;
  private currentViewer: DG.DockNode | null = null;
  
  constructor(df: DG.DataFrame, ligands: DG.Column, target: string, poses: number) {
    this.df = df;
    this.ligands = ligands;
    this.target = target;
    this.poses = poses;
  }
  
  async run() {
    const { posesColumn, confidenceColumn, virtualPosesColumn } = await this.calculatePoses();
    this.posesColumn = posesColumn;
    this.virtualPosesColumn = virtualPosesColumn;
    this.subscribeToCurrentCellChanged();
  }
  
  private async calculatePoses() {
    const posesColumn = this.createColumn(DG.TYPE.STRING, CONSTANTS.POSES_COLUMN_NAME, this.ligands.length);
    const confidenceColumn = this.createColumn(DG.TYPE.INT, CONSTANTS.CONFIDENCE_COLUMN_NAME, this.ligands.length);
    const virtualPosesColumn = this.createColumn(DG.TYPE.OBJECT, CONSTANTS.VIRTUAL_POSES_COLUMN_NAME, this.ligands.length);
  
    for (let i = 0; i < posesColumn.length; ++i) {
      const posesJson = await this.getPosesJson(this.ligands.get(i));
      const pdbHelper: IPdbHelper = await getPdbHelper();
      const { bestPose, confidence } = this.findBestPose(posesJson);
  
      posesColumn.set(i, await pdbHelper.molToPdb(bestPose));
      posesColumn.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
      //needed in
      posesColumn.setTag('docking.role', 'ligand');
      confidenceColumn.set(i, confidence);
      posesJson.receptor = this.target;
      virtualPosesColumn.set(i, posesJson);
    }
  
    this.df.columns.add(posesColumn);
    this.df.columns.add(confidenceColumn);
    this.df.columns.add(virtualPosesColumn);
    await grok.data.detectSemanticTypes(this.df);
  
    return { posesColumn, confidenceColumn, virtualPosesColumn };
  }
  
  public async getPosesJson(ligand: string): Promise<PosesJson> {
    const encodedPoses = await grok.functions.call('Bionemo:diffdock', {
      protein: this.target,
      ligand: ligand,
      num_poses: this.poses,
    });
    return JSON.parse(new TextDecoder().decode(encodedPoses.data));
  }
  
  public createColumn(type: DG.ColumnType, name: string, length: number): DG.Column {
    const unusedName = this.df.columns.getUnusedName(name);
    return DG.Column.fromType(type, unusedName, length);
  }
  
  private findBestPose(poses: PosesJson) {
    const positions = poses.position_confidence;
    const minValue = Math.min(...positions);
    const idx = positions.findIndex((value) => value === minValue);
    return {
      bestPose: poses.ligand_positions[idx],
      confidence: minValue,
    };
  }
  
  private createPoseItems(posesJson: PosesJson): string[] {
    return posesJson.ligand_positions.map((_, index) =>
      `Pose ${index + 1} - Score: ${posesJson.position_confidence[index]}`
    );
  }
  
  public async createCombinedControl(currentRow: number) {
    const molstarViewer = await this.createMolstarViewer(currentRow);
    molstarViewer.root.style.height = '100%';
    const items = this.createPoseItems(this.virtualPosesColumn.get(currentRow) as PosesJson);
    
    const molstarProps = molstarViewer.getProperties();
    const ligandValueProp = molstarProps.find((p: DG.Property) => p.name === 'ligandValue');
    const ligandColProp = molstarProps.find((p: DG.Property) => p.name === 'ligandColumnName');
  
    const comboPopup = ui.comboPopup(items[0], items, (item: string) => {
      comboPopup.getElementsByTagName('span')[0].textContent = item;
      const idx = parseInt(item.split(' ')[1]) - 1;
      molstarViewer.apply({ 'ligandValue': (this.virtualPosesColumn.get(currentRow) as PosesJson).ligand_positions[idx] });
      molstarViewer.onPropertyChanged(ligandValueProp!);
      molstarViewer.onPropertyChanged(ligandColProp!);
    });
  
    comboPopup.classList.add('bionemo-combo-popup');
    return ui.divV([comboPopup, molstarViewer.root]);
  }
  
  private async createMolstarViewer(currentRow: number) {
    return await this.posesColumn.dataFrame.plot.fromType('Biostructure', {
      pdb: this.target,
      ligandColumnName: this.ligands.name,
      ligandValue: (this.virtualPosesColumn.get(currentRow) as PosesJson).ligand_positions[currentRow],
      zoom: true,
    });
  }

  private async handleCurrentCellChanged() {
    const currentCell = this.df.currentCell;
    if (currentCell && currentCell.column === this.posesColumn) {
      const currentRow = this.df.currentRowIdx;
      const combinedControl = await this.createCombinedControl(currentRow);
      const view = grok.shell.getTableView(this.df.name);
  
      if (this.currentViewer)
        view.dockManager.close(this.currentViewer);

      this.currentViewer = view.dockManager.dock(ui.divV([combinedControl]), DG.DOCK_TYPE.RIGHT);
    }
  }

  private subscribeToCurrentCellChanged() {
    this.df.onCurrentCellChanged.subscribe(() => this.handleCurrentCellChanged());
  }
}