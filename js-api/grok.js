import * as _chem from './src/chem';
import * as _ml from './src/ml';
import {Dapi} from "./src/dapi";
import {Functions} from "./src/functions";
import {Events} from "./src/events";
import {Settings, Shell} from "./src/shell";
import {Data} from "./src/data";


export const functions = new Functions();
export const events = new Events();
export const dapi = new Dapi();
export const shell = new Shell();
export const settings = new Settings();
export const data = new Data();
export const chem = _chem;
export const ml = _ml;
