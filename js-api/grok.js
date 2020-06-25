import * as _chem from './src/chem';
import * as _ml from './src/ml';
import {Dapi} from "./src/dapi";
import {Functions} from "./src/functions";
import {Events} from "./src/events";
import {Settings, Shell} from "./src/shell";
import {Data} from "./src/data";


export let functions = new Functions();
export let events = new Events();
export let dapi = new Dapi();
export let shell = new Shell();
export let settings = new Settings();
export let data = new Data();
export let chem = _chem;
export let ml = _ml;
