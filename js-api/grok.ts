import {Dapi} from './src/dapi';
import {Functions} from './src/functions';
import {Events} from './src/events';
import {Settings, Shell} from './src/shell';
import {Data} from './src/data';
import {Logger} from './src/logger';
import {UserSettingsStorage} from "./src/user_settings_storage";


export const functions = new Functions();
export const events = new Events();
export const dapi = new Dapi();
export const shell = new Shell();
export const settings = new Settings();
export const data = new Data();
export const userSettings = new UserSettingsStorage();
export * from './src/chem';
export * from './src/ml';
export * from './src/decorators/functions';

export const log = Logger.getStatic();