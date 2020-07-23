import * as _chem from './src/chem';
import * as _ml from './src/ml';
import {Dapi} from "./src/dapi";
import {Functions} from "./src/functions";
import {Events} from "./src/events";
import {Settings, Shell} from "./src/shell";
import {Data} from "./src/data";

export let functions: Functions;
export let events: Events;
export let dapi: Dapi;
export let shell: Shell;
export let settings: Settings;
export let data: Data;
export let chem: typeof _chem;
export let ml: typeof _ml;