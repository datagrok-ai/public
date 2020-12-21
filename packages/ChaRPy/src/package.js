/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import mapR from "./mapR.json";
import mapPy from "./mapPy.json";

export let _package = new DG.Package();

// Recursive substitution function #1
// Takes in a multi-level getOptions() object and replaces
// all the user inputs with R/Python alternatives from mapR/Py.json
// Input: target (type: Object), parameters received from getOptions()
//        source (type: Object), R/Python substitutes
// Output: target (type: Object), a modified parameters object with all values substituted
function assignOnlyIntersection(target, source) {
  Object.keys(target).forEach(key => {
    (target[key] == null) && delete target[key]
    if (target[key] in source) {
      target[key] = source[target[key]];
    } else if (typeof target[key] === 'object') {
      assignOnlyIntersection(target[key], source)
    }
  })
}

// Recursive substitution function #2
// Takes in a generalized R/Python viewer code and slices it according to user selected options
// Input: groupByList (type: list), an empty list that will be filled only with column names
//        colsFilter
//        stRing (type: string), a generalized R/Python viewer code obtained from map.json
//        optionsObj (type: Object), getOptions() output
//        map (type: Object), a predefined mapping of all user selected parameters to R/Python alternatives
// Output: stRing (type: string), complete R/Python plot script
//         groupByList (type: list), contains column names
//         colsFilter
function dynamicReplace(groupByList, colsFilter, stRing, optionsObj, map) {
  Object.keys(optionsObj).forEach(key => {
    (optionsObj[key] == null) && delete optionsObj[key];

    if (key.includes('ColumnName')) {
      if (key !== 'valueColumnName') {
        groupByList.push(optionsObj[key]);
      }
      if (typeof optionsObj[key] === 'object') {
        Object.keys(optionsObj[key]).forEach(k => {
          colsFilter.push(optionsObj[key][k]);
        })
      } else {
        colsFilter.push(optionsObj[key]);
      }
    }

    if (typeof optionsObj[key] === 'object' &&
      Object.keys(optionsObj[key]).every(elem => elem !== '0')) {
      stRing = dynamicReplace(groupByList, colsFilter, stRing, optionsObj[key], map);
    } else {
      stRing = stRing.split("!(" + key + ")").join(map[key]);
      stRing = stRing.split("!(" + key + ")").join(optionsObj[key]);
    }
  })
  return stRing;
}

// Table preprocessing function
// Created a new truncated dataframe
// Input: colsFilter (type: list), list of columns to keep
//        table (type: dataframe), original dataframe
// Output: t (type: dataframe), new truncated dataframe
function tablePreProcess(colsFilter, table) {
  let l = [];
  for (let j = 0; j < colsFilter.length; j++) {
    if (colsFilter[j] !== '') {
      l.push(table.columns.byName(colsFilter[j]));
    }
  }
  let t = DG.DataFrame.fromColumns(l);
  return t;
}

// Top-level string substitution function that implements both recursive functions
// Carries out a sequence of substitutions to produce a finalized script for the selected environment
// Input: optionsObj (type: Object), getOptions() output
//        map (type: Object), a predefined mapping of all user selected parameters to R/Python alternatives
// Output: stRing (type: string), complete R plot script
async function strReplace(table, optionsObj, map, language) {

  let string;
  let paramsMap;

  // extract the generalized viewer code and the corresponding parameter mapping
  if (optionsObj.type === 'Trellis plot') {

    if (language === 'R') {
      string = map.plotScripts[optionsObj.look.innerViewerLook['#type'].toLowerCase().replace('look', '')] + " + facet_grid(!(xColumnNames)!(yColumnNames))";
    } else if (language === 'Python') {
      string = map.plotScripts[optionsObj.look.innerViewerLook['#type'].toLowerCase().replace('look', '')];
    }
    paramsMap = map.additionalOps[optionsObj.look.innerViewerLook['#type'].toLowerCase().replace('look', '')];

  } else {
    string = map.plotScripts[optionsObj.type.toLowerCase().replace(' ', '')];
    paramsMap = map.additionalOps[optionsObj.type.toLowerCase().replace(' ', '')];
  }

  // adjust the getOptions() output by replacing misc grok codes with R analogues
  assignOnlyIntersection(optionsObj, map.miscCodes);

  // adjust the generalized viewer code by substituting in the values from
  // mapR.json and getOptions() output
  let groupByList = [];
  let colsFilter = [];
  string = dynamicReplace(groupByList, colsFilter, string, optionsObj, paramsMap);
  string = string.replace("!(groupByList)", groupByList);
  for (let id = 0; id < table.columns.length; id++) {
    if (table.columns.byIndex(id).type === 'datetime') {
      string = string.replace("!(dateConvert)",
        map.dateConvert.split("!(dateConvert)").join(table.columns.names()[id]));
    }
  }
  string = string.replace(/!\([^)]*\) */g, "");

  let viewerLeft = DG.Viewer.fromType(optionsObj.type, table, optionsObj.look);
  let viewerRight = DG.Viewer.fromType('Scripting Viewer',
    tablePreProcess([...new Set(colsFilter)], table),
    {script: map.header + string + map.tail});

  // create a container for viewers
  let block =
    $(ui.splitV([
      ui.textArea(string),
      ui.splitH([
        viewerLeft,
        viewerRight])
    ])).css('flex-grow', '1');

  // create and show a dialogue window containing the generated code string,
  // original viewer and the scripting viewer output
  let dialog = ui.dialog('OUTPUT SCRIPT').add(block[0]);
  dialog.onClose.subscribe((_) => {
    viewerLeft.dataFrame = new DG.DataFrame();
    viewerRight.dataFrame = new DG.DataFrame();
  });
  dialog.showModal(true);
}

// Creates menu buttons that executes viewer to code conversion
grok.events.onContextMenu.subscribe((args) => {
  if (args.args.context instanceof DG.Viewer) {
    let menu = args.args.menu.group('To Script');

    menu.item('to R', async () => {
      let options = JSON.parse(args.args.context.getOptions());
      await strReplace(args.args.context.table, options, mapR, 'R')
    });
    menu.item('to Python', async () => {
      let options = JSON.parse(args.args.context.getOptions());
      await strReplace(args.args.context.table, options, mapPy, 'Python')
    });
  }
});

//name: exportFunc
//tags: autostart
export function toScriptInit() {
}