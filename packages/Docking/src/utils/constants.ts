import * as DG from 'datagrok-api/dg';
import { DockingPackage } from '../package-utils';

export const _package = new DockingPackage();
export const TARGET_PATH = 'System:AppData/Docking/targets';
export const CACHED_RESULTS: DG.LruCache<string, DG.DataFrame> = new DG.LruCache<string, DG.DataFrame>();
export let BINDING_ENERGY_COL = 'binding energy';
export let BINDING_ENERGY_COL_UNUSED = '';
export let POSE_COL = 'pose';
export let POSE_COL_UNUSED = '';
export const ERROR_COL_NAME = 'error';
export const ERROR_MESSAGE = 'Error log';

export function setPose(value: string) {
    POSE_COL_UNUSED = value;
}

export function setAffinity(value: string) {
    BINDING_ENERGY_COL_UNUSED = value;
}

export const PROPERTY_DESCRIPTIONS: {[colName: string]: string} = {
    'intermolecular (1)': 'Final Intermolecular Energy',
    'electrostatic': 'Electrostatic Energy',
    'ligand fixed': 'Moving Ligand-Fixed Receptor',
    'ligand moving': 'Moving Ligand-Moving Receptor',
    'total internal (2)': 'Final Total Internal Energy',
    'torsional free (3)': 'Torsional Free Energy',
    'unbound systems (4)': 'Unbound System\s Energy' 
}