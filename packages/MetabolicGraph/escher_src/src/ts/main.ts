/* eslint-disable camelcase */
/**
* @license
*
* Escher
* Author: Zachary King
*
* The MIT License (MIT)
*
* This software is Copyright © 2015 The Regents of the University of
* California. All Rights Reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*/

// ESCHER_VERSION provided by webpack plugin or main-node
/* global ESCHER_VERSION */

// --JSX

import underscore from 'underscore';
import preact from 'preact';
import * as baconjs from 'baconjs';
import mousetrap from 'mousetrap';
import vkbeautify from 'vkbeautify';
import {select as d3_select, selection as d3_selection} from 'd3-selection';
import {json as d3_json} from 'd3-request';

export const version = '1.0.0';

// @ts-ignore
export {default as Builder} from '../Builder';
export {EscherMap} from './escherMap';
export {default as Behavior} from './Behavior';
export {default as KeyManager} from './KeyManager';
export {default as DataMenu} from './DataMenu';
export {default as UndoStack} from './UndoStack';
export {default as CobraModel} from './CobraModel';
export * as utils from './utils';
export {default as SearchIndex} from './SearchIndex';
export {default as Settings} from './Settings';
export * as dataStyles from './dataStyles';
export {default as ZoomContainer} from './ZoomContainer';


export const libs = {
  _: underscore,
  underscore,
  preact,
  baconjs,
  mousetrap,
  vkbeautify,
  d3_selection,
  d3_select,
  d3_json
};
