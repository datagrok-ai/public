/* eslint-disable camelcase */
/* global Blob, XMLSerializer, Image, btoa */

import vkbeautify from 'vkbeautify';
import _ from 'underscore';
import {json as d3_json} from 'd3-request';
import {text as d3_text} from 'd3-request';
import {csvParseRows as d3_csvParseRows} from 'd3-dsv';
import {Selection as d3_selection} from 'd3-selection';
import {Coord, D3Selection, D3SVGSelection} from './types';

// const vkbeautify = require('vkbeautify');
// const _ = require('underscore');
// const d3_json = require('d3-request').json;
// const d3_text = require('d3-request').text;
// const d3_csvParseRows = require('d3-dsv').csvParseRows;

try {
  // @ts-ignore
  if (!window.saveAs)
    // @ts-ignore
    window.saveAs = require('file-saver').saveAs;
} catch (e) {
  console.warn('Not a browser, so FileSaver.js not available.');
}

/**
 * Check if Blob is available, and alert if it is not.
 */
export function _check_filesaver() {
  try {
    const isFileSaverSupported = !!new Blob();
  } catch (e) {
    alert('Blob not supported');
  }
}

export function set_options<T extends Object>(options: Partial<T>, defaults: T, must_be_float?: { [key: string]: any }): T {
  if (options === undefined || options == null)
    return defaults as T;

  const i = -1;
  const out: { [key: string]: any } = {};
  for (const key in defaults) {
    const has_key = ((key in options) &&
                   (options[key] !== null) &&
                   (options[key] !== undefined));
    let val: any = (has_key ? options[key] : defaults[key]);
    if (must_be_float && key in must_be_float) {
      val = parseFloat(val);
      if (isNaN(val)) {
        if (has_key) {
          console.warn('Bad float for option ' + key);
          val = parseFloat(defaults[key] as string);
          if (isNaN(val)) {
            console.warn('Bad float for default ' + key);
            val = null;
          }
        } else {
          console.warn('Bad float for default ' + key);
          val = null;
        }
      }
    }
    out[key] = val;
  }
  return out as T;
}


export function remove_child_nodes(selection: D3Selection) {
  /** Removes all child nodes from a d3 selection

   */
  const node = selection.node();
  while (node?.hasChildNodes())
    node.removeChild(node.lastChild as Node);
}


export function load_css(css_path: string, callback: (css: string) => void) {
  let css = '';
  if (css_path) {
    d3_text(css_path, function(error, text) {
      if (error)
        console.warn(error);

      css = text;
      callback(css);
    });
  }
  return false;
}


export function _ends_with(str: string, suffix: string) {
  return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

/**
 * Load a file.
 * @param {} t - this context for callback. Should be an object.
 * @param {} files_to_load - A filename to load. Must be JSON or CSS.
 * @param {} callback - Function to run after the file is loaded. Takes the
 * arguments error and data.
 * @param {} value - If the value is specified, just assign it and do not
 * execute the ajax query.
 */
export function load_the_file(t: any, file: string, callback: (error: any, d: any) => void, value?: any) {
  if (value) {
    if (file) console.warn('File ' + file + ' overridden by value.');
    callback.call(t, null, value);
    return;
  }
  if (!file) {
    callback.call(t, 'No filename', null);
    return;
  }
  if (_ends_with(file, 'json'))
    d3_json(file, function(e, d) { callback.call(t, e, d); });
  else if (_ends_with(file, 'css'))
    d3_text(file, function(e, d) { callback.call(t, e, d); });
  else
    callback.call(t, 'Unrecognized file type', null);

  return;
}


export function load_files(t: any, files_to_load: { file: string, callback: (error: any, d: any) => void, value: any }[], final_callback: (error: any, d: any) => void) {
  /** Load multiple files asynchronously by calling utils.load_the_file.

      t: this context for callback. Should be an object.

      files_to_load: A list of objects with the attributes:

      { file: a_filename.json, callback: a_callback_fn }

      File must be JSON or CSS.

      final_callback: Function that runs after all files have loaded.

  */
  if (files_to_load.length === 0) final_callback.call(t, null, null);
  let i = -1;
  let remaining = files_to_load.length;
  while (++i < files_to_load.length) {
    load_the_file(
      t,
      files_to_load[i].file,
      function(e: any, d: any) {
        // @ts-ignore
        (this as any).call(t, e, d);
        if (!--remaining) final_callback.call(t, null, null);
      }.bind(files_to_load[i].callback),
      files_to_load[i].value
    );
  }
}


/**
 * Create a constructor that returns a new object with our without the 'new'
 * keyword.
 *
 * Adapted from Hubert Kauker (MIT Licensed), John Resig (MIT Licensed).
 * http://stackoverflow.com/questions/7892884/simple-class-instantiation
 */
export function make_class() {
  let is_internal: boolean;
  const constructor = function(args: any) {
    // @ts-ignore
    // eslint-disable-next-line no-invalid-this
    if (this instanceof constructor) {
      // @ts-ignore
      // eslint-disable-next-line no-invalid-this
      if (typeof this.init === 'function')
        // @ts-ignore
        // eslint-disable-next-line prefer-rest-params, no-invalid-this
        this.init.apply(this, ...(is_internal ? args : arguments as any));
    } else {
      is_internal = true;
      // @ts-ignore
      // eslint-disable-next-line prefer-rest-params
      const instance = new constructor(arguments as any);
      is_internal = false;
      return instance;
    }
  };
  return constructor;
}

/**
 * Return a class that can be instantiated without the new keyword.
 * @param {Class} AClass - Any ES6 class.
 */
export function class_with_optional_new(AClass: any) {
  return new Proxy(AClass, {
    apply(Target, thisArg, args) {
      return new Target(...args);
    }
  });
}

export function setup_defs(svg: D3Selection, style: string) {
  // add stylesheet
  svg.select('defs').remove();
  const defs = svg.append('defs');
  // make sure the defs is the first node
  const node = defs.node()!;
  node.parentNode!.insertBefore(node, node.parentNode!.firstChild);
  defs.append('style')
    .attr('type', 'text/css')
    .text(style);
  return defs;
}

/**
 * Run through the d3 data binding steps for an object. Also checks to make sure
 * none of the values in the *object* are undefined, and ignores those.
 *
 * The create_function, update_function, and exit_function CAN modify the input
 * data object.
 *
 * @param {} container_sel - A d3 selection containing all objects.
 *
 * @param {} parent_node_selector - A selector string for a subselection of
 * container_sel.
 *
 * @param {} children_selector - A selector string for each DOM element to bind.
 *
 * @param {} object - An object to bind to the selection.
 *
 * @param {} id_key - The key that will be used to store object IDs in the bound
 * data points.
 *
 * @param {} create_function - A function for enter selection. Create function
 * must return a selection of the new nodes.
 *
 * @param {} update_function - A function for update selection.
 *
 * @param {} exit_function - A function for exit selection.
 */
export function draw_an_object(container_sel: D3Selection, parent_node_selector: string, children_selector: string,
  object: any, id_key: string, create_function: any, update_function?: any,
  exit_function?: any) {
  const draw_object: { [key: string]: any } = {};

  for (const id in object) {
    if (object[id] === undefined)
      console.warn('Undefined value for id ' + id + ' in object. Ignoring.');
    else
      draw_object[id] = object[id];
  }

  const sel = container_sel.select(parent_node_selector)
    .selectAll(children_selector)
    .data(make_array_ref(draw_object, id_key),
      function(d: any) { return d[id_key]; });

  // enter: generate and place reaction
  const update_sel = create_function ?
    create_function(sel.enter()).merge(sel) :
    sel;

  // update: update when necessary
  if (update_function)
    update_sel.call(update_function);


  // exit
  if (exit_function)
    sel.exit().call(exit_function);
}

/**
 * Run through the d3 data binding steps for an object that is nested within
 * another element with D3 data.
 *
 * The create_function, update_function, and exit_function CAN modify the input
 * data object.
 *
 * @param {} container_sel - A d3 selection containing all objects.
 *
 * @param {} children_selector - A selector string for each DOM element to bind.
 *
 * @param {} object_data_key - A key for the parent object containing data for
 * the new selection.
 *
 * @param {} id_key - The key that will be used to store object IDs in the bound
 * data points.
 *
 * @param {} create_function - A function for enter selection. Create function
 * must return a selection of the new nodes.
 *
 * @param {} update_function - A function for update selection.
 *
 * @param {} exit_function - A function for exit selection.
 */
export function draw_a_nested_object(container_sel: D3Selection, children_selector: string, object_data_key: string,
  id_key: string, create_function: any, update_function: any,
  exit_function: any) {
  const sel = container_sel.selectAll(children_selector)
    .data(function(d: any) { return make_array_ref(d[object_data_key], id_key); },
      function(d: any) { return d[id_key]; });

  // enter: generate and place reaction
  const update_sel = (create_function ?
    create_function(sel.enter()).merge(sel) :
    sel);

  // update: update when necessary
  if (update_function)
    update_sel.call(update_function);


  // exit
  if (exit_function)
    sel.exit().call(exit_function);
}

export function make_array(obj: any, id_key: string) { // is this super slow?
  const array = [];
  for (const key in obj) {
    // copy object
    const it = clone(obj[key]);
    // add key as 'id'
    it[id_key] = key;
    // add object to array
    array.push(it);
  }
  return array;
}

export function make_array_ref(obj: any, id_key: string) {
  /** Turn the object into an array, but only by reference. Faster than
      make_array.

  */
  const array = [];
  for (const key in obj) {
    // copy object
    const it = obj[key];
    // add key as 'id'
    it[id_key] = key;
    // add object to array
    array.push(it);
  }
  return array;
}

export function compare_arrays(a1: any, a2: any) {
  /** Compares two simple (not-nested) arrays.

   */
  if (!a1 || !a2) return false;
  if (a1.length != a2.length) return false;
  for (let i = 0, l = a1.length; i < l; i++) {
    if (a1[i] != a2[i]) {
      // Warning - two different object instances will never be equal: {x:20} != {x:20}
      return false;
    }
  }
  return true;
}

/**
 * Convert an array of objects to an object with all keys and values
 * that are arrays of the same length as arr. Fills in spaces with null.
 *
 *  For example, [ { a: 1 }, { b: 2 }] becomes { a: [1, null], b: [null, 2] }.
 */
export function arrayToObject(arr: any) {
  // new object
  const obj: { [key: string]: any } = {};
  // for each element of the array
  for (let i = 0, l = arr.length; i < l; i++) {
    const column = arr[i];
    const keys = Object.keys(column);
    for (let k = 0, nk = keys.length; k < nk; k++) {
      const id = keys[k];
      if (!(id in obj)) {
        const n = [];
        // fill spaces with null
        for (let j = 0; j < l; j++)
          n[j] = null;

        n[i] = column[id];
        obj[id] = n;
      } else
        obj[id][i] = column[id];
    }
  }
  return obj;
}

/**
 * Deep copy for array and object types. All other types are returned by
 * reference.
 * @param {T<Object|Array|*>} obj - The object to copy.
 * @return {T} The copied object.
 */
export function clone<T = any>(obj: T): T {
  if (_.isArray(obj))
    return _.map(obj, function(t: any) { return clone(t); }) as T;
  else if (_.isObject(obj))
    return _.mapObject(obj, function(t: any, k: any) { return clone(t); }) as T;
  else
    return obj;
}

export function extend<T extends Object = any >(obj1: T, obj2: Partial<T>, overwrite: boolean = false) {
  /** Extends obj1 with keys/values from obj2. Performs the extension
      cautiously, and does not override attributes, unless the overwrite
      argument is true.

      Arguments
      ---------

      obj1: Object to extend

      obj2: Object with which to extend.

      overwrite: (Optional, Default false) Overwrite attributes in obj1.

  */

  if (overwrite === undefined)
    overwrite = false;

  for (const attrname in obj2) {
    if (!(attrname in obj1) || overwrite) // UNIT TEST This
      obj1[attrname] = obj2[attrname] as any;
    else
      throw new Error('Attribute ' + attrname + ' already in object.');
  }
}

export function uniqueConcat(arrays: any[]) {
  const newArray: any[] = [];
  arrays.forEach((a: any[]) => {
    a.forEach((x: any) => {
      if (newArray.indexOf(x) < 0)
        newArray.push(x);
    });
  });
  return newArray;
}

/**
 * Return unique values in array of strings.
 *
 * http://stackoverflow.com/questions/1960473/unique-values-in-an-array
 */
export function unique_strings_array(arr: any[]) {
  const a = [];
  for (let i = 0, l = arr.length; i < l; i++) {
    if (a.indexOf(arr[i]) === -1)
      a.push(arr[i]);
  }
  return a;
}

/**
 * Returns a function, that, as long as it continues to be invoked, will not be
 * triggered. The function will be called after it stops being called for N
 * milliseconds. If "immediate" is passed, trigger the function on the leading
 * edge, instead of the trailing.
 */
export function debounce(func: Function, wait: number, immediate: boolean) {
  let timeout: any;
  return function() {
    // @ts-ignore
    // eslint-disable-next-line no-invalid-this
    const context = this as any;
    // @ts-ignore
    // eslint-disable-next-line prefer-rest-params
    const args = arguments as any;
    const later = function() {
      timeout = null;
      if (!immediate) func.apply(context, args);
    };
    const callNow = immediate && !timeout;
    clearTimeout(timeout);
    timeout = setTimeout(later, wait);
    if (callNow) func.apply(context, args);
  };
}

/**
 * Return a copy of the object with just the given ids.
 * @param {} obj - An object
 * @param {} ids - An array of id strings
 */
export function object_slice_for_ids(obj: any, ids: string[]) {
  const subset: { [key: string]: any } = {};
  let i = -1;
  while (++i < ids.length)
    subset[ids[i]] = clone(obj[ids[i]]);

  if (ids.length !== Object.keys(subset).length)
    console.warn('did not find correct reaction subset');

  return subset;
}

/**
 * Return a reference of the object with just the given ids. Faster than
 * object_slice_for_ids.
 * @param {} obj - An object.
 * @param {} ids - An array of id strings.
 */
export function object_slice_for_ids_ref(obj: any, ids: string[]) {
  const subset: { [key: string]: any } = {};
  let i = -1;
  while (++i < ids.length)
    subset[ids[i]] = obj[ids[i]];

  if (ids.length !== Object.keys(subset).length)
    console.warn('did not find correct reaction subset');

  return subset;
}

export function c_plus_c(coords1: Coord | null | undefined, coords2: Coord | null | undefined) {
  if (coords1 == null || coords2 == null ||
      coords1 === undefined || coords2 === undefined)
    return null;

  return {
    x: coords1.x + coords2.x,
    y: coords1.y + coords2.y,
  };
}

export function c_minus_c(coords1: Coord | null | undefined, coords2: Coord | null | undefined) {
  if (coords1 == null || coords2 == null ||
      coords1 === undefined || coords2 === undefined)
    return null;

  return {
    x: coords1.x - coords2.x,
    y: coords1.y - coords2.y,
  };
}

export function c_times_scalar(coords: Coord, scalar: number) {
  return {
    x: coords.x * scalar,
    y: coords.y * scalar,
  };
}

/**
 * Download JSON file in a blob.
 */
export function download_json(json: any, name: string) {
  // Alert if blob isn't going to work
  _check_filesaver();

  const j = JSON.stringify(json);
  const blob = new Blob([j], {type: 'application/json'});
  window.saveAs?.(blob, name + '.json');
}

/**
 * Try to load the file as JSON.
 * @param {} f - The file path
 * @param {} callback - A callback function that accepts arguments: error, data.
 * @param {} pre_fn (optional) - A function to call before loading the data.
 * @param {} failure_fn (optional) - A function to call if the load fails or is
 * aborted.
*/
export function load_json(f: string, callback: (error: any, d: any) => void, pre_fn: () => void, failure_fn: () => void) {
  // Check for the various File API support
  if (!(window.File && window.FileReader && window.FileList && window.Blob))
    callback('The File APIs are not fully supported in this browser.', null);


  const reader = new window.FileReader();
  // Closure to capture the file information.
  reader.onload = function(event: any) {
    const result = event.target.result;
    let data;
    // Try JSON
    try {
      data = JSON.parse(result);
    } catch (e) {
      // If it failed, return the error
      callback(e, null);
      return;
    }
    // If successful, return the data
    callback(null, data);
  };
  if (pre_fn !== undefined && pre_fn !== null)
    try { pre_fn(); } catch (e) { console.warn(e); }

  reader.onabort = function(event) {
    try { failure_fn(); } catch (e) { console.warn(e); }
  };
  reader.onerror = function(event) {
    try { failure_fn(); } catch (e) { console.warn(e); }
  };
  // Read in the image file as a data URL
  reader.readAsText(new Blob([f]));
}

/**
 * Try to load the file as JSON or CSV (JSON first).
 * @param {String} f - The file path
 * @param {Function}  csv_converter - A function to convert the CSV output to equivalent JSON.
 * @param {Function} callback - A callback function that accepts arguments: error, data.
 * @param {} pre_fn (optional) - A function to call before loading the data.
 * @param {} failure_fn (optional) - A function to call if the load fails or is
 * aborted.
 * @param {} debug_event (optional) - An event, with a string at
 * event.target.result, to load as though it was the contents of a loaded file.
 */
export function load_json_or_csv(f: string | File | Blob, csv_converter: (data: any) => any, callback: (error: any, d: any) => void, pre_fn?: () => void, failure_fn?: () => void,
  debug_event?: any) {
  // Capture the file information.
  const onload_function = function(event: any) {
    const result = event.target.result;
    let data;
    let errors;
    // try JSON
    try {
      data = JSON.parse(result);
    } catch (e) {
      errors = 'JSON error: ' + e;

      // try csv
      try {
        data = csv_converter(d3_csvParseRows(result));
      } catch (e) {
        // if both failed, return the errors
        callback(errors + '\nCSV error: ' + e, null);
        return;
      }
    }
    // if successful, return the data
    callback(null, data);
  };
  if (debug_event !== undefined && debug_event !== null) {
    console.warn('Debugging load_json_or_csv');
    return onload_function(debug_event);
  }

  // Check for the various File API support.
  if (!(window.File && window.FileReader && window.FileList && window.Blob))
    callback('The File APIs are not fully supported in this browser.', null);
  const reader = new window.FileReader();

  if (pre_fn !== undefined && pre_fn !== null)
    try { pre_fn(); } catch (e) { console.warn(e); }

  reader.onabort = function(event) {
    try { failure_fn(); } catch (e) { console.warn(e); }
  };
  reader.onerror = function(event) {
    try { failure_fn(); } catch (e) { console.warn(e); }
  };
  // Read in the image file as a data URL.
  reader.onload = onload_function;
  reader.readAsText(new Blob([f]));
}

/**
 * Download an svg file using FileSaver.js.
 * @param {String} name - The filename (without extension)
 * @param {D3 Selection} svg_sel - The d3 selection for the SVG element
 * @param {Boolean} do_beautify - If true, then beautify the SVG output
 */
export function downloadSvg(name: string, svg_sel: D3Selection, do_beautify: boolean) {
  // Alert if blob isn't going to work
  _check_filesaver();

  // Make the xml string
  let xml = (new XMLSerializer()).serializeToString(svg_sel.node()!);
  if (do_beautify) xml = vkbeautify.xml(xml);
  xml = ('<?xml version="1.0" encoding="utf-8"?>\n' +
         '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n' +
         ' "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n' +
         xml);

  // Save
  const blob = new Blob([xml], {type: 'image/svg+xml'});
  window.saveAs?.(blob, name + '.svg');
}

/**
 * Download a png file using FileSaver.js.
 * @param {String} name - The filename (without extension).
 * @param {D3 Selection} svg_sel - The d3 selection for the SVG element.
 */
export function downloadPng(name: string, svg_sel: D3SVGSelection) {
  // Alert if blob isn't going to work
  _check_filesaver();

  // Make the xml string
  let xml = new XMLSerializer().serializeToString(svg_sel.node()!);
  xml = ('<?xml version="1.0" encoding="utf-8"?>\n' +
         '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"\n' +
         ' "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n' +
         xml);

  // Canvas to hold the image
  const canvas = document.createElement('canvas');
  const context = canvas.getContext('2d')!;

  // Get SVG size
  const svg_size = svg_sel.node()!.getBBox();
  const svg_width = svg_size.width + svg_size.x;
  const svg_height = svg_size.height + svg_size.y;

  // Canvas size = SVG size. Constrained to 10000px for very large SVGs
  if (svg_width < 10000 && svg_height < 10000) {
    canvas.width = svg_width;
    canvas.height = svg_height;
  } else {
    if (canvas.width > canvas.height) {
      canvas.width = 10000;
      canvas.height = 10000 * (svg_height / svg_width);
    } else {
      canvas.width = 10000 * (svg_width / svg_height);
      canvas.height = 10000;
    }
  }

  // Image element appended with data
  const base_image = new Image();
  base_image.src = 'data:image/svg+xml;base64,' + btoa(xml);

  base_image.onload = function() {
    // Draw image to canvas with white background
    context.fillStyle = '#FFF';
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(base_image, 0, 0, canvas.width, canvas.height);

    // Save image
    canvas.toBlob(function(blob: Blob | null) {
      blob && window.saveAs?.(blob, name + '.png');
    });
  };
}

export function rotate_coords_recursive(coords_array: Coord[], angle: number, center: Coord) {
  return coords_array.map(function(c) {
    return rotate_coords(c, angle, center);
  });
}

/**
 * Calculates displacement { x: dx, y: dy } based on rotating point c around
 * center with angle.
 */
export function rotate_coords(c: Coord, angle: number, center: Coord) {
  const dx = Math.cos(-angle) * (c.x - center.x) +
      Math.sin(-angle) * (c.y - center.y) +
    center.x - c.x;
  const dy = - Math.sin(-angle) * (c.x - center.x) +
    Math.cos(-angle) * (c.y - center.y) +
    center.y - c.y;
  return {x: dx, y: dy};
}

/**
 * Get the angle between coordinates
 * @param {Object} coords - Array of 2 coordinate objects { x: 1, y: 1 }
 * @return {Number} angle between 0 and 2PI.
 */
export function get_angle(coords: Coord[]) {
  const denominator = coords[1].x - coords[0].x;
  const numerator = coords[1].y - coords[0].y;
  if (denominator === 0 && numerator >= 0)
    return Math.PI / 2;
  else if (denominator === 0 && numerator < 0)
    return 3 * Math.PI / 2;
  else if (denominator >= 0 && numerator >= 0)
    return Math.atan(numerator / denominator);
  else if (denominator >= 0)
    return (Math.atan(numerator / denominator) + 2 * Math.PI);
  else
    return (Math.atan(numerator / denominator) + Math.PI);
}

export function to_degrees(radians: number) {
  return radians * 180 / Math.PI;
}

/**
 * Force to domain -PI to PI
 */
export function angleNorm(radians: number) {
  if (radians < -Math.PI)
    return radians + Math.floor((radians - Math.PI) / (-2 * Math.PI)) * 2 * Math.PI;
  else if (radians > Math.PI)
    return radians - Math.floor((radians + Math.PI) / (2 * Math.PI)) * 2 * Math.PI;
  else
    return radians;
}

export function to_radians(degrees: number) {
  return Math.PI / 180 * degrees;
}

/**
 * Convert to radians, and force to domain -PI to PI
 */
export function to_radians_norm(degrees: number) {
  const radians = to_radians(degrees);
  return angleNorm(radians);
}

export function angle_for_event(displacement: Coord, point: Coord, center: Coord) {
  const gamma = Math.atan2((point.x - center.x), (center.y - point.y));
  const beta = Math.atan2((point.x - center.x + displacement.x),
    (center.y - point.y - displacement.y));
  const angle = beta - gamma;
  return angle;
}

export function distance(start: Coord, end: Coord) {
  return Math.sqrt(Math.pow(end.y - start.y, 2) + Math.pow(end.x - start.x, 2));
}

/**
 * Report an error if any of the arguments are undefined. Call by passing in
 * "arguments" from any function and an array of argument names.
 */
export function check_undefined(args: any[], names: string[]) {
  names.forEach(function(name, i) {
    if (args[i] === undefined)
      console.error(`Argument is undefined: ${names[i]}`);
  });
}

export function compartmentalize(bigg_id: string, compartment_id: string) {
  return `${bigg_id}_${compartment_id}`;
}

/**
 * Returns an array of [bigg_id, compartment id]. Matches compartment ids with
 * length 1 or 2. Return [ id, null ] if no match is found.
 */
export function decompartmentalize(id: string) {
  const reg = /(.*)_([a-z0-9]{1,2})$/;
  const result = reg.exec(id);
  return result !== null ? result.slice(1, 3) : [id, null];
}

export function mean(array: number[]) {
  const sum = array.reduce(function(a, b) { return a + b; });
  const avg = sum / array.length;
  return avg;
}

export function median(array: number[]) {
  array.sort(function(a, b) { return a - b; });
  const half = Math.floor(array.length / 2);
  if (array.length % 2 == 1)
    return array[half];
  else
    return (array[half - 1] + array[half]) / 2.0;
}

export function quartiles(array: number[]) {
  array.sort(function(a, b) { return a - b; });
  const half = Math.floor(array.length / 2);
  if (array.length === 1) {
    return [
      array[0],
      array[0],
      array[0],
    ];
  } else if (array.length % 2 === 1) {
    return [
      median(array.slice(0, half)),
      array[half],
      median(array.slice(half + 1)),
    ];
  } else {
    return [
      median(array.slice(0, half)),
      (array[half - 1] + array[half]) / 2.0,
      median(array.slice(half)),
    ];
  }
}

/**
 * Generate random characters
 *
 * Thanks to @csharptest.net
 * http://stackoverflow.com/questions/1349404/generate-a-string-of-5-random-characters-in-javascript
 */
export function random_characters(num: number) {
  let text = '';
  const possible = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
  for (let i = 0; i < num; i++)
    text += possible.charAt(Math.floor(Math.random() * possible.length));

  return text;
}

export function generate_map_id() {
  return random_characters(12);
}

/**
 * Check that the selection has the given parent tag.
 * @param {D3 Selection|DOM Node} el - A D3 Selection or DOM Node to check.
 * @param {String} tag - A tag name (case insensitive).
 */
export function check_for_parent_tag(el: D3Selection | Element, tag: string) {
  // make sure it is a node
  if (!(el instanceof Element))
    el = el.node()!;

  while (el.parentNode != null) {
    el = el.parentNode as Element;
    if (el.tagName == undefined)
      continue;

    if (el.tagName.toLowerCase() === tag.toLowerCase())
      return true;
  }
  return false;
}

/**
 * Convert model or map name to url.
 * @param {String} name - The short name, e.g. e_coli.iJO1366.central_metabolism.
 * @param {String} download_url (optional) - The url to prepend.
 */
export function name_to_url(name: string, download_url: string | null | undefined) {
  if (download_url !== undefined && download_url !== null) {
    // strip download_url
    download_url = download_url.replace(/^\/|\/$/g, '');
    name = [download_url, name].join('/');
  }
  // strip final path
  return name.replace(/^\/|\/$/g, '') + '.json';
}

/**
 * Get the document for the node
 */
export function get_document(node: Element) {
  return node.ownerDocument;
}

/**
 * Get the window for the node
 */
export function get_window(node: Element) {
  return get_document(node).defaultView;
}

/**
 * Get translation and rotation values for a transform string. This used to be
 * in d3, but since v4, I just adapted a solution from SO:
 *
 * http://stackoverflow.com/questions/38224875/replacing-d3-transform-in-d3-v4
 *
 * To get skew and scale out, go back to that example.
 *
 * TODO rename function without "catch"
 *
 * @param {String} transform_attr - A transform string.
 */
export function d3_transform_catch(transform_attr: string) {
  if (transform_attr.indexOf('skew') !== -1 ||
      transform_attr.indexOf('matrix') !== -1)
    throw new Error('d3_transform_catch does not work with skew or matrix');


  const translate_res = (/translate\s*\(\s*([0-9.-]+)\s*,\s*([0-9.-]+)\s*\)/
    .exec(transform_attr));
  const tn = _.isNull(translate_res);
  const tx = tn ? 0.0 : Number(translate_res[1]);
  const ty = tn ? 0.0 : Number(translate_res[2]);

  const rotate_res = (/rotate\s*\(\s*([0-9.-]+)\s*\)/
    .exec(transform_attr));
  const rn = _.isNull(rotate_res);
  const r = rn ? 0.0 : Number(rotate_res[1]);

  const scale_res = (/scale\s*\(\s*([0-9.-]+)\s*\)/
    .exec(transform_attr));
  const sn = _.isNull(scale_res);
  const s = sn ? 0.0 : Number(scale_res[1]);

  return {translate: [tx, ty], rotate: r, scale: s,};

  // // Create a dummy g for calculation purposes only. This will new be appended
  // // to the DOM and will be discarded once this function returns.
  // var g = document.createElementNS('http://www.w3.org/2000/svg', 'g')

  // // Set the transform attribute to the provided string value.
  // g.setAttributeNS(null, 'transform', transform_attr)

  // // Consolidate the SVGTransformList containing all Try to a single
  // // SVGTransform of type SVG_TRANSFORM_MATRIX and get its SVGMatrix.

  // var matrix = g.transform.baseVal.consolidate().matrix

  // // Below calculations are taken and adapted from the private func
  // // transform/decompose.js of D3's module d3-interpolate.
  // var a = matrix.a
  // var b = matrix.b
  // var c = matrix.c
  // var d = matrix.d
  // var e = matrix.e
  // var f = matrix.f
  // var scaleX = Math.sqrt(a * a + b * b)

  // if (scaleX) {
  //   a /= scaleX
  //   b /= scaleX
  // }

  // if (a * d < b * c) {
  //   a = -a
  //   b = -b
  // }

  // return {
  //   translate: [ e, f ],
  //   rotate: Math.atan2(b, a) * Math.PI / 180,
  // }
}

/**
 * Look for name in the user agent string.
 */
// function check_browser (name) {
//   var browser = function() {
//     // Thanks to
//     // http://stackoverflow.com/questions/2400935/browser-detection-in-javascript
//     var ua = navigator.userAgent
//     var M = ua.match(/(opera|chrome|safari|firefox|msie|trident(?=\/))\/?\s*(\d+)/i) || []
//     var tem
//     if (/trident/i.test(M[1])) {
//       tem = /\brv[ :]+(\d+)/g.exec(ua) || []
//       return 'IE '+ (tem[1] || '')
//     }
//     if (M[1] === 'Chrome') {
//       tem = ua.match(/\b(OPR|Edge)\/(\d+)/)
//       if (tem != null) return tem.slice(1).join(' ').replace('OPR', 'Opera')
//     }
//     M = M[2] ? [ M[1], M[2] ]: [ navigator.appName, navigator.appVersion, '-?' ]
//     if ((tem = ua.match(/version\/(\d+)/i)) !== null) {
//       M.splice(1, 1, tem[1])
//     }
//     return M.join(' ')
//   }

//   try {
//     // navigator.userAgent is deprecated, so don't count on it
//     return browser().toLowerCase().indexOf(name) > -1
//   } catch (e) {
//     return false
//   }
// }
