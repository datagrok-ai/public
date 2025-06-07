import { PlateQuery, PropertyCondition } from "../plates-crud";

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from "rxjs/operators";
import {NumericMatcher, StringInListMatcher} from '../numeric_matcher';
import {PlateProperty, plateProperties, getPlateUniquePropertyValues, getWellUniquePropertyValues, queryPlates, queryWells, wellProperties} from "../plates-crud";


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

function getPlatesSearchForm(): DG.InputForm {
  return getSearchForm(plateProperties, getPlateUniquePropertyValues);
}

function getWellsSearchForm(): DG.InputForm {
  return getSearchForm(wellProperties, getWellUniquePropertyValues);
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
  const platesForm = getPlatesSearchForm();
  const wellsForm = getWellsSearchForm();

  const refresh = () => {
    const query: PlateQuery = {
      plateMatchers: searchFormToMatchers(platesForm),
      wellMatchers: searchFormToMatchers(wellsForm)
    };

    search(query).then((df) => {
      view.dataFrame = df;
      onResults(view.grid);
    });
  }

  platesForm.onInputChanged.pipe(debounceTime(500)).subscribe((_) => refresh());
  wellsForm.onInputChanged.pipe(debounceTime(500)).subscribe((_) => refresh());

  //@ts-ignore
  view.dockManager.dock(ui.divH([platesForm.root, wellsForm.root]), DG.DOCK_TYPE.TOP);

  refresh();
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