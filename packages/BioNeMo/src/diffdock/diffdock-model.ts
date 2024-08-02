import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {IPdbHelper, getPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import '../../css/bionemo.css';

interface PosesJson {
  ligand_positions: string[];
  position_confidence: number[];
  receptor: string;
}
  
export const CONSTANTS = {
  POSES_COLUMN_NAME: 'poses',
  CONFIDENCE_COLUMN_NAME: 'confidence',
  VIRTUAL_POSES_COLUMN_NAME: '~posesJson'
};
  
export class DiffDockModel {
  private df: DG.DataFrame;
  private ligands: DG.Column;
  private target: string;
  private poses: number;
  public posesColumn!: DG.Column;
  public virtualPosesColumn!: DG.Column;
  private currentViewer: DG.DockNode | null = null;
  private bestId: number | null = null;
  
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
    const confidenceColumn = this.createColumn(DG.TYPE.FLOAT, CONSTANTS.CONFIDENCE_COLUMN_NAME, this.ligands.length);
    const virtualPosesColumn = this.createColumn(DG.TYPE.STRING, CONSTANTS.VIRTUAL_POSES_COLUMN_NAME, this.ligands.length);
    const grid = grok.shell.getTableView(this.df.name).grid;

    for (let i = 0; i < posesColumn.length; ++i) {
      const posesJson = await this.getPosesJson(this.ligands.get(i));
      const pdbHelper: IPdbHelper = await getPdbHelper();
      const { bestPose, confidence } = this.findBestPose(posesJson);
  
      posesColumn.set(i, await pdbHelper.molToPdb(bestPose));
      posesColumn.setTag(DG.TAGS.SEMTYPE, DG.SEMTYPE.MOLECULE3D);
      //needed in order not to open custom molstar
      posesColumn.setTag('docking.role', 'ligand');

      confidenceColumn.set(i, confidence);
      confidenceColumn.setTag(DG.TAGS.FORMAT, '0.00');
      confidenceColumn.meta.colors.setLinear([DG.Color.green, DG.Color.red]);
      grid.sort([confidenceColumn]);

      posesJson.receptor = this.target;
      virtualPosesColumn.set(i, JSON.stringify(posesJson));
    }
  
    this.df.columns.add(posesColumn);
    this.df.columns.add(confidenceColumn);
    this.df.columns.add(virtualPosesColumn);
    await grok.data.detectSemanticTypes(this.df);

    grid.onCellRender.subscribe((args: any) => {
      grid.setOptions({ 'rowHeight': 100 });
      grid.col(CONSTANTS.POSES_COLUMN_NAME)!.width = 100;
      grid.col(CONSTANTS.CONFIDENCE_COLUMN_NAME)!.width = 100 + 50;
      grid.col(CONSTANTS.CONFIDENCE_COLUMN_NAME)!.isTextColorCoded = true;
      grid.invalidate();
    });
  
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
    this.bestId = idx;
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
  
  public async createCombinedControl(currentRow: number, fullHeight: boolean = true): Promise<HTMLDivElement> {
    const molstarViewer = await this.createMolstarViewer(currentRow);
    if (fullHeight)
      molstarViewer.root.style.height = '100%';
    const items = this.createPoseItems(JSON.parse(this.virtualPosesColumn.get(currentRow)) as PosesJson);
    
    const molstarProps = molstarViewer.getProperties();
    const ligandValueProp = molstarProps.find((p: DG.Property) => p.name === 'ligandValue');
    const ligandColProp = molstarProps.find((p: DG.Property) => p.name === 'ligandColumnName');
  
    const comboPopup = ui.comboPopup(items[this.bestId ?? 0], items, (item: string) => {
      comboPopup.getElementsByTagName('span')[0].textContent = item;
      const idx = parseInt(item.split(' ')[1]) - 1;
      molstarViewer.apply({ 'ligandValue': (JSON.parse(this.virtualPosesColumn.get(currentRow)) as PosesJson).ligand_positions[idx] });
      molstarViewer.onPropertyChanged(ligandValueProp!);
      molstarViewer.onPropertyChanged(ligandColProp!);
    });
  
    comboPopup.classList.add('bionemo-combo-popup');
    const container = ui.divV([comboPopup, molstarViewer.root]);
    container.style.height = '100%';
    return container;
  }
  
  private async createMolstarViewer(currentRow: number): Promise<DG.Widget> {
    return await this.ligands.dataFrame.plot.fromType('Biostructure', {
      pdb: this.target,
      ligandColumnName: this.ligands.name,
      ligandValue: (JSON.parse(this.virtualPosesColumn.get(currentRow)) as PosesJson).ligand_positions[this.bestId ?? 0],
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

  public subscribeToCurrentCellChanged() {
    this.df.onCurrentCellChanged.subscribe(() => this.handleCurrentCellChanged());
  }
}