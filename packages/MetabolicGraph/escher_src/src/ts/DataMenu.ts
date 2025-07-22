/* eslint-disable camelcase */
/** DataMenu */
import * as utils from './utils';
import {json as d3_json} from 'd3-request';
import {D3Selection, SettingsType} from './types';

export default function DataMenu(options: {selection: D3Selection, data: any, datafiles: any, getdatafiles: any, update_callback: any, target: any}) {
  const o: {selection: D3Selection, data: any, datafiles: any, getdatafiles: any, update_callback: any, target: any} = utils.set_options(options, {
    selection: null,
    getdatafiles: null,
    datafiles: null,
    update_callback: null,
    target: null,
    data: null});

  if (o.selection == null)
    throw new Error('No selection provided for DataMenu');

  // setup dropdown menu
  // Append menu if it doesn't exist
  let menu = o.selection.select('.data-menu');
  if (menu.empty()) {
    menu = o.selection.append('div')
      .attr('class', 'data-menu');
  }
  const select_sel = menu.append('form')
    .append('select').attr('class', 'dropdown-menu');

  if (o.getdatafiles) {
    if (o.datafiles)
      console.warn('DataMenu: getdatafiles option overrides datafiles');

    d3_json(o.getdatafiles, function(error, d) {
      // returns json object:  { data: [file0, file1, ...] }
      if (error)
        return console.warn(error);
      else
        load_with_files(o.target, d.data, select_sel, o.update_callback, o.selection);

      return null;
    });
  } else if (o.datafiles)
    load_with_files(o.target, o.datafiles, select_sel, o.update_callback, o.selection);
  else
    console.warn('DataMenu: No datafiles given');


  return {update: update};

  // definitions
  function load_with_files(t: any, files: any, select_sel: any, update_callback: any, selection: any) {
    //when it changes
    select_sel.node().addEventListener('change', function() {
      // @ts-ignore
      // eslint-disable-next-line no-invalid-this
      load_datafile(t, this.value, selection, update_callback);
    }, false);

    const file = files[0];

    update(files, select_sel);
    load_datafile(t, file, selection, update_callback);
  };
  function load_datafile(t: any, this_file: any, selection: any, callback: any) {
    utils.load_the_file(t, this_file, function(error: any, data: any) {
      if (error) {
        selection.append('error loading');
        o.data = null;
        return console.warn(error);
      } else {
        o.data = data;
        if (callback)
          callback(data);
      }
    });
  };

  function update(list: any, select_sel: any) {
    // update select element with d3 selection /select_sel/ to have options
    // given by /list/
    // TODO remove paths from file list
    select_sel.selectAll('.menu-option')
      .data(list)
      .enter()
      .append('option')
      .attr('value', function(d: any) { return d; } )
      .text(function(d: any) { return d; } );
    // TODO set value to default_filename_index
    select_sel.node().focus();
  };

  function get_data() {
    return o.data;
  };
};
