import {Shell} from './src/shell';
import {Logger} from './src/logger';

/** Visual shell (projects, panels, views) */
export const shell = new Shell();

export const log = Logger.getStatic();

export * from './src/chem';
