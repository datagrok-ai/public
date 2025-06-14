import { PlateQuery, PlateTemplate, PropertyCondition } from "../plates-crud";

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from "rxjs/operators";
import {NumericMatcher, StringInListMatcher} from '../numeric_matcher';
import {PlateProperty, plateProperties, getPlateUniquePropertyValues, getWellUniquePropertyValues, queryPlates, queryWells, wellProperties, plateTemplates} from "../plates-crud";


type PropInput = DG.InputBase & { prop: PlateProperty };


function getSearchForm(properties: PlateProperty[], getUniqueValues: (prop: any) => string[]): DG.InputForm {
  const inputs: PropInput[] = [];
  for (const prop of properties) {
    if (prop.value_type == DG.TYPE.STRING) {
      const items = getUniqueValues(prop);
      const input = ui.input.multiChoice(prop.name, { items: items}) as PropInput;
      input.prop = prop;
      inputs.push(input);
    }
    else if (prop.value_type == DG.TYPE.FLOAT) {
      const input = ui.input.string(prop.name) as PropInput;
      input.prop = prop;
      input.addValidator(s => NumericMatcher.parse(s) ? null : 'Invalid numerical criteria. Example: ">10"');
      inputs.push(input);
    }
  }

  return DG.InputForm.forInputs(inputs);
}

function getPlatesSearchForm(filteredProperties?: PlateProperty[]): DG.InputForm {
  return getSearchForm(filteredProperties ?? plateProperties, getPlateUniquePropertyValues);
}

function getWellsSearchForm(filteredProperties?: PlateProperty[]): DG.InputForm {
  return getSearchForm(filteredProperties ?? wellProperties, getWellUniquePropertyValues);
}

function searchFormToMatchers(form: DG.InputForm): PropertyCondition[] {
  const matchers: PropertyCondition[] = [];

  for (const input of form.inputs as PropInput[]) {
    if (input.inputType == DG.InputType.MultiChoice && input.value.length > 0) {
      matchers.push({
        property: input.prop,
        matcher: new StringInListMatcher(input.value)
      });
    }
    else if (input.inputType == DG.InputType.Text && input.value !== '' && NumericMatcher.parse(input.value)) {
      matchers.push({
        property: input.prop,
        matcher: NumericMatcher.parse(input.value)!
      });
    }
  }

  return matchers;
}

function getSearchView(search: (query: PlateQuery) => Promise<DG.DataFrame>, onResults: (grid: DG.Grid) => void): DG.View {
  const dummy = DG.DataFrame.create(5, 'Search plates');
  const view = DG.TableView.create(dummy);
  let platesForm = getPlatesSearchForm();
  let wellsForm = getWellsSearchForm();
  const platesFormHost = ui.div();
  const wellsFormHost = ui.div();
  let plateTemplate = plateTemplates[0];

  const refreshUI = () => {
    ui.empty(platesFormHost);
    ui.empty(wellsFormHost);
    platesForm = getPlatesSearchForm(plateTemplate?.plateProperties?.map(p => p as PlateProperty));
    wellsForm = getWellsSearchForm(plateTemplate?.wellProperties?.map(p => p as PlateProperty));
    platesForm.onInputChanged.pipe(debounceTime(500)).subscribe((_) => refreshResults());
    wellsForm.onInputChanged.pipe(debounceTime(500)).subscribe((_) => refreshResults());
    platesFormHost.appendChild(platesForm.root);
    wellsFormHost.appendChild(wellsForm.root);
  }

  const setTemplate = (template: PlateTemplate) => {
    plateTemplate = template;
    refreshUI();
  }

  const plateTemplateSelector = ui.input.choice('Template', {
    nullable: true,
    items: plateTemplates.map(pt => pt.name),
    value: plateTemplates[0].name,
    onValueChanged: (v) => setTemplate(plateTemplates.find(pt => pt.name === v)!)
  });


  const refreshResults = () => {
    const query: PlateQuery = {
      plateMatchers: searchFormToMatchers(platesForm),
      wellMatchers: searchFormToMatchers(wellsForm)
    };

    search(query).then((df) => {
      view.dataFrame = df;
      onResults(view.grid);
    });
  }

  const searchHost = ui.divV([
    plateTemplateSelector.root,
    ui.divH([
      ui.divV([ui.h2('Plates'), platesFormHost], {style: {flexGrow: '1'}}), 
      ui.divV([ui.h2('Wells'), wellsFormHost], {style: {flexGrow: '1'}})
    ], {style: {height: '100%', width: '100%'}}),
  ], {classes: 'ui-panel'});
  view.dockManager.dock(searchHost, DG.DOCK_TYPE.TOP);

  //refreshUI();
  setTemplate(plateTemplates[0]);
  refreshResults();
  return view;
}

export function searchPlatesView(): DG.View {
  return getSearchView(queryPlates, (grid) => {
    grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
        .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`
  });
}

export function searchWellsView(): DG.View {
  return getSearchView(queryWells, (grid) => {});
}