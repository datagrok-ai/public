import {Dapi} from './src/dapi';
import {Functions} from './src/functions';
import {Events} from './src/events';
import {Settings, Shell} from './src/shell';
import {Data} from './src/data';
import {Logger} from './src/logger';
import {UserSettingsStorage} from "./src/user_settings_storage";
import {AI} from "./src/ai";

/** Function-related APIs (finding, calling, registering) */
export const functions = new Functions();

/** Global events (listening, firing) */
export const events = new Events();

/** Server API (entities, queries, functions, credentials, jobs, users, projects, notebooks, models, packages, layouts, views, tables, groups, scripts, environments) */
export const dapi = new Dapi();

/** Visual shell (projects, panels, views, viewers, current objects, etc) */
export const shell = new Shell();

/** Settings API (settings, shell, dapi, userSettings) */
export const settings = new Settings();

/** Creating, loading, querying, manipulating, joining tables. */
export const data = new Data();

/** User settings API (userSettings, userSettingsStorage) */
export const userSettings = new UserSettingsStorage();

/** AI API */
export const ai = AI;

export * from './src/chem';
export * from './src/ml';
export * from './src/decorators/functions';

/** Generated REST client (@ExportApi server routers); call {@link dapi2Init}(grok.dapi.root, grok.dapi.token) first. */
export {dapi2, init as dapi2Init} from './src/dapi.g';

export const log = Logger.getStatic();
