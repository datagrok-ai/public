
//import $ from 'cash-dom';
import * as dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
dayjs.extend(utc);

import * as wu from 'wu';

import * as _const from './src/const';
import * as events from './src/events';
import * as dapi from './src/dapi';
import * as dataframe from './src/dataframe';
import * as entities from './src/entities';
import * as ddd_api_g from './src/api/ddt.api.g';
import * as grok_shared_api_g from './src/api/grok_shared.api.g';
import * as shell from './src/shell';
import * as functions from './src/functions';
import * as wrappers from './src/wrappers_impl';
import {time, timeAsync, Utils, HtmlUtils, LruCache} from './src/utils';
import * as sticky_meta from './src/sticky_meta';
import * as data from './src/data';
import * as helpers from './src/helpers';
import * as logger from './src/logger';

import * as chem from './src/chem';
import * as ml from './src/ml';
import * as utils from './src/utils';

import * as grok from './grok';

export let DG = {
    ...utils,
    ...ml,
    ...chem,
    ...logger,
    ...helpers,
    ...data,
    ...sticky_meta,
    ..._const,
    ...dapi,
    ...events,
    ...dataframe,
    ...entities,
    ...ddd_api_g,
    ...grok_shared_api_g,
    ...shell,
    ...wrappers,
    ...functions,
    time, timeAsync, Utils, HtmlUtils, LruCache,
};

export * as grok from './grok';

(DG as any).DG = DG;
(globalThis as any).grok = grok;
(globalThis as any).DG = DG;

const api: any = (typeof window !== 'undefined' ? window : global.window) as any;
api.dayjs = dayjs;
api.wu = wu;
api.grok = grok;
(globalThis as any).wu = wu;
(globalThis as any).dayjs = dayjs;

