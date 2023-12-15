import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {
  AutoDockRunResult, GridSize, IAutoDockService
} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {BiostructureData, BiostructureDataJson} from '@datagrok-libraries/bio/src/pdb/types';
import {getPdbHelper, IPdbHelper} from '@datagrok-libraries/bio/src/pdb/pdb-helper';
import {DockerContainerStatus, awaitStatus} from '@datagrok-libraries/bio/src/utils/docker';


import {_package, getAutoDockService} from '../package';
import {DockerContainer} from 'datagrok-api/dg';
import {delay} from '@datagrok-libraries/utils/src/test';

namespace Forms {
  export type run = {
    /** PDB string of the receptor target structure */ receptor: string,
    /** PDB string of the ligand structure */ ligand: string, /* TODO: BiostructureDataJson support */
    /** */ poseCount?: number,
  }
  export type runRes = {
    /** pdbqt */ poses: string;
  }
}

export function buildDefaultGridConfig(receptorName: string, gridSize: GridSize): string {
  const res: string = `
npts ${gridSize.x} ${gridSize.y} ${gridSize.y}      # num.grid points in xyz
gridfld ${receptorName}.maps.fld                    # grid_data_file
spacing 0.375                                       # spacing(A)
receptor_types A C Fe N NA OA SA                    # receptor atom types
ligand_types A C N F NA OA HD SA                    # ligand atom types
receptor ${receptorName}.pdbqt                      # macromolecule 
gridcenter auto                                     # xyz-coordinates or auto
smooth 0.5                                          # store minimum energy w/in rad(A)
map ${receptorName}.A.map                           # atom-specific affinity map
map ${receptorName}.C.map                           # atom-specific affinity map
map ${receptorName}.N.map                           # atom-specific affinity map
map ${receptorName}.F.map                           # atom-specific affinity map
map ${receptorName}.NA.map                          # atom-specific affinity map
map ${receptorName}.OA.map                          # atom-specific affinity map
map ${receptorName}.HD.map                          # atom-specific affinity map
map ${receptorName}.SA.map                          # atom-specific affinity map
elecmap ${receptorName}.e.map                       # electrostatic potential map
dsolvmap ${receptorName}.d.map                      # desolvation potential map
dielectric -0.1465                                  # <0, AD4 distance-dep.diel;>0, constant
`.split('\n').filter((s) => !!s).join('\n');
  return res;
}

export class AutoDockService implements IAutoDockService {
  private readonly dcName: string;
  private dc!: DG.DockerContainer;
  private ph: IPdbHelper;

  constructor() {
    this.dcName = `${_package.name.toLowerCase()}-autodock`;
  }

  async startDockerContainer(timeout: number = 30000): Promise<void> {
    // TODO: Use the new dockerContainers API
    const res = await grok.dapi.docker.dockerContainers.run(this.dc.id /*, true */);
    let end: boolean = false;
    for (let i = 0; i < timeout / 200; ++i) {
      this.dc = await grok.dapi.docker.dockerContainers.find(this.dc.id);
      switch (this.dc.status) {
        case 'stopped': {
          await grok.dapi.docker.dockerContainers.run(this.dc.id);
          break;
        }
        case 'pending change':
        case 'changing': {
          // skip to wait
          break;
        }
        case 'checking':
        case 'started': {
          end = true;
          break;
        }
        case 'error': {
          throw new Error('Docker container error state.');
        }
      }
      if (end) break;
      await delay(200);
    }
    if (!end) throw new Error('Docker container run timeout.');
    this.dc = await grok.dapi.docker.dockerContainers.find(this.dc.id);
  }

  async init(): Promise<void> {
    [this.dc, this.ph] = await Promise.all([
      // returns docker container with status not actual soon
      grok.dapi.docker.dockerContainers.filter(this.dcName).first(),
      getPdbHelper()
    ]);
    if (!this.dc)
      throw new Error(`AutoDock docker container not found '${this.dcName}'.`);
  }

  get ready(): boolean {
    return this.dc.status === 'started' || this.dc.status === 'checking';
  }

  async awaitStatus(targetStatus: DockerContainerStatus, timeout: number = 30000): Promise<void> {
    return awaitStatus(this.dc.id, targetStatus, timeout, _package.logger);
  }

  async run(
    receptor: string, ligand: BiostructureData, npts: GridSize, poseCount?: number,
    poseColName?: string, debug: boolean = false
  ): Promise<AutoDockRunResult> {
    if (ligand.binary || ligand.ext !== 'pdb')
      throw new Error(`Unsupported ligand ext '${ligand.ext}' or binary, must be 'pdb' string.`);
    const ligandPdb = ligand.data as string;

    const form: Forms.run = {
      receptor: receptor,
      ligand: ligandPdb,
      poseCount: poseCount,
    };
    const params: RequestInit = {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(form),
    };

    const path = `/dock?x=${npts.x}&y=${npts.y}&z=${npts.z}&debug=${debug}`;
    const adResStr = (await grok.dapi.docker.dockerContainers.request(this.dc.id, path, params))!;
    // TODO: Use the new dockerContainers API
    // const adResponse: Response = await grok.dapi.docker.dockerContainers.fetchProxy(this.dc.id, path, params);
    // if (adResponse.status !== 200) {
    //   const errMsg = (await adResponse.json())['datagrok-error'];
    //   throw new Error(errMsg);
    // }
    //const adRes: Forms.runRes = (await adResponse.json()) as Forms.runRes;
    const adRes: Forms.runRes = JSON.parse(adResStr) as Forms.runRes;

    // const modelList: string[] = wu(adRes.poses.matchAll(/MODEL.*?ENDMDL/gs/* lazy, not greedy */))
    //   .map((ma) => ma[0]).toArray();
    // const posesCol = DG.Column.fromStrings(poseColName ?? 'pdbqt_model', modelList);
    // // posesCol.semType = DG.SEMTYPE.MOLECULE3D;
    // // posesCol.setTag(DG.TAGS.UNITS, 'pdbqt');
    // const posesDf = DG.DataFrame.fromColumns([posesCol]);

    const posesDf: DG.DataFrame = this.ph.parsePdbqt(adRes.poses, poseColName);

    const res: AutoDockRunResult = {
      posesDf: posesDf,
    };
    return res;
  }

  static async getSvc(): Promise<IAutoDockService> {
    const svc: AutoDockService = new AutoDockService();
    await svc.init();
    return svc;
  }
}

export async function _runAutodock(
  receptor: DG.FileInfo, ligand: DG.FileInfo, x: number, y: number, z: number
): Promise<AutoDockRunResult | null> {
  const svc: IAutoDockService = new AutoDockService();
  if (!svc.ready) {
    grok.shell.warning('Autodock container not started yet.');
    return null;
  }

  const receptorStr = await receptor.readAsString();
  const ligandData: BiostructureData = {binary: true, ext: ligand.extension, data: await ligand.readAsBytes()};
  return await svc.run(receptorStr, ligandData, new GridSize(x, y, z));
}

export async function _runAutodock2(molCol: DG.Column<string>, receptorPdb: string): Promise<void> {
  // const receptorPdb: string = await receptorFi.readAsString();

  let resDf: DG.DataFrame | undefined = undefined;

  const adSvc = await getAutoDockService();
  await adSvc.awaitStatus('started', 30000);
  const ph = await getPdbHelper();
  for (let lRowI = 0; lRowI < molCol.length; ++lRowI) {
    const t1 = window.performance.now();
    try {
      const ligandMol = molCol.get(lRowI);
      const ligandPdb = await ph.molToPdb(ligandMol!);
      const ligandData: BiostructureData = {binary: false, data: ligandPdb, ext: 'pdb'};

      const adRes = await adSvc.run(
        receptorPdb, ligandData, new GridSize(40, 40, 40), 10);

      if (resDf === undefined)
        resDf = adRes.posesDf.clone(DG.BitSet.create(adRes.posesDf.rowCount, (_i) => false));

      resDf!.append(adRes.posesDf, true);
    } finally {
      const t2 = window.performance.now();
      _package.logger.debug('_runAutodock2(), ' + `ligand: ${lRowI}, ` + `ET: ${t2 - t2} ms, `);
    }
  }

  if (resDf !== undefined)
    DG.Utils.download('models.csv', resDf.toCsv());
  else
    window.alert('Empty result');
}

