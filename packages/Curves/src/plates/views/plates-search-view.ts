/* eslint-disable max-len */
import {getPlateUniquePropertyValues, getWellUniquePropertyValues, PlateQuery, PlateTemplate, PropertyCondition} from '../plates-crud';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {debounceTime} from 'rxjs/operators';
import {NumericMatcher, StringInListMatcher} from '../numeric_matcher';
import * as rxjs from 'rxjs';
import {
  PlateProperty,
  getUniquePropertyValues, // We'll use this single function for both plates and wells
  queryPlates,
  queryWells,
  allProperties, // Import the new unified cache
  plateTemplates
} from '../plates-crud';


type PropInput = DG.InputBase & { prop: PlateProperty };


function getSearchForm(
  allProperties: PlateProperty[],
  templateProperties: PlateProperty[],
  // CHANGED: The function signature is simplified, as the underlying unique value logic is now unified.
  getUniqueValues: (prop: any) => string[]
): { root: HTMLElement, templateForm: DG.InputForm, otherForm: DG.InputForm } {
  const createInput = (prop: PlateProperty): PropInput | null => {
  // Check if choices exist (inherited from IProperty)
    if (prop.choices && prop.choices.length > 0) {
    // If choices is a string (from DB), parse it; otherwise use as-is
      let choicesList: string[];
      if (typeof prop.choices === 'string') {
        try {
          choicesList = JSON.parse(prop.choices as any);
        } catch {
        // If parsing fails, treat as a single choice
          choicesList = [prop.choices as any];
        }
      } else {
        choicesList = prop.choices;
      }

      const input = ui.input.multiChoice(prop.name, {items: choicesList}) as PropInput;
      input.prop = prop;
      return input;
    } else if (prop.type == DG.TYPE.STRING) {
      const input = ui.input.multiChoice(prop.name, {items: getUniqueValues(prop)}) as PropInput;
      input.prop = prop;
      return input;
    } else if (prop.type == DG.TYPE.FLOAT || prop.type == DG.TYPE.INT) {
      const input = ui.input.string(prop.name) as PropInput;
      input.prop = prop;
      input.addValidator((s) => NumericMatcher.parse(s) ? null : 'Invalid numerical criteria. Example: ">10"');
      return input;
    }
    return null;
  }; ; ;

  const templatePropNames = new Set(templateProperties.map((p) => p.name));

  const templateInputs = templateProperties.map(createInput).filter((i) => i) as PropInput[];
  // FIXED: Filter from the full list of properties, not just globals.
  const otherInputs = allProperties.filter((p) => !templatePropNames.has(p.name))
    .map(createInput).filter((i) => i) as PropInput[];

  const templateForm = DG.InputForm.forInputs(templateInputs);
  const otherForm = DG.InputForm.forInputs(otherInputs);

  const templateContainer = ui.divV([
    ui.h3('Template Properties', {style: {marginTop: '0px', marginBottom: '2px'}}),
    templateForm.root
  ], {
    style: {
      border: '1px solid var(--grey-2)',
      borderRadius: '5px',
      padding: '8px',
      marginBottom: '10px'
    }
  });

  const root = ui.divV([
    ...(templateInputs.length > 0 ? [templateContainer] : []),
    ...(otherInputs.length > 0 ? [otherForm.root] : [])
  ]);

  return {root, templateForm, otherForm};
}


function searchFormToMatchers(form: DG.InputForm): PropertyCondition[] {
  const matchers: PropertyCondition[] = [];
  if (!form) return matchers;

  for (const input of form.inputs as PropInput[]) {
    if (input.inputType == DG.InputType.MultiChoice && input.value.length > 0) {
      matchers.push({
        property: input.prop,
        matcher: new StringInListMatcher(input.value)
      });
    } else if (input.inputType == DG.InputType.Text && input.value !== '' && NumericMatcher.parse(input.value)) {
      matchers.push({
        property: input.prop,
        matcher: NumericMatcher.parse(input.value)!
      });
    }
  }
  return matchers;
}


function getSearchView(
  search: (query: PlateQuery) => Promise<DG.DataFrame>,
  onResults: (grid: DG.Grid) => void
): DG.View {
  const mainView = DG.View.create();
  const resultsView = DG.TableView.create(DG.DataFrame.create(0), false);

  let platesTemplateForm: DG.InputForm; let platesOtherForm: DG.InputForm;
  let wellsTemplateForm: DG.InputForm; let wellsOtherForm: DG.InputForm;

  const platesFormHost = ui.div();
  const wellsFormHost = ui.div();
  let plateTemplate = plateTemplates[0];

  const plateTemplateSelector = ui.input.choice('Template', {
    items: plateTemplates.map((pt) => pt.name),
    value: plateTemplates[0].name,
    onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)
  });

  const globalSearchInput = ui.input.bool('Search all properties', {
    value: false,
    onValueChanged: () => {
      plateTemplateSelector.enabled = !globalSearchInput.value;
      refreshUI();
    }
  });
  globalSearchInput.setTooltip('Search using all properties in the database, not just those defined in the template.');

  const refreshResults = () => {
    const query: PlateQuery = {
      plateMatchers: [...searchFormToMatchers(platesTemplateForm), ...searchFormToMatchers(platesOtherForm)],
      wellMatchers: [...searchFormToMatchers(wellsTemplateForm), ...searchFormToMatchers(wellsOtherForm)]
    };
    search(query).then((df) => {
      resultsView.dataFrame = df;
      onResults(resultsView.grid);
    });
  };

  const refreshUI = () => {
    ui.empty(platesFormHost);
    ui.empty(wellsFormHost);

    const currentTemplatePlateProps = plateTemplate?.plateProperties?.map((p) => p as PlateProperty) ?? [];
    const currentTemplateWellProps = plateTemplate?.wellProperties?.map((p) => p as PlateProperty) ?? [];

    // REWRITTEN: This is the main fix.
    // Instead of using the old variables, we now filter the unified 'allProperties' cache by scope.
    // In plates-search-view.ts, update the refreshUI function:

    const allPlateProps = globalSearchInput.value ?
      Array.from(new Map(
        allProperties
          .filter((p) => p.scope === 'plate')
          .map((p) => [p.name, p]) // Use name as key to deduplicate
      ).values()) :
      currentTemplatePlateProps;

    const allWellProps = globalSearchInput.value ?
      Array.from(new Map(
        allProperties
          .filter((p) => p.scope === 'well')
          .map((p) => [p.name, p]) // Use name as key to deduplicate
      ).values()) :
      currentTemplateWellProps;

    // FIXED: Pass the simplified 'getUniquePropertyValues' function for both calls.
    // NEW, CORRECTED CODE
    const platesFormsResult = getSearchForm(allPlateProps, currentTemplatePlateProps, getPlateUniquePropertyValues);
    const wellsFormsResult = getSearchForm(allWellProps, currentTemplateWellProps, getWellUniquePropertyValues);

    platesTemplateForm = platesFormsResult.templateForm;
    platesOtherForm = platesFormsResult.otherForm;
    wellsTemplateForm = wellsFormsResult.templateForm;
    wellsOtherForm = wellsFormsResult.otherForm;

    [platesTemplateForm.root, platesOtherForm.root, wellsTemplateForm.root, wellsOtherForm.root].forEach((r) => r.style.width = '100%');

    const combinedSub = DG.debounce(
      rxjs.merge(
        platesTemplateForm.onInputChanged,
        platesOtherForm.onInputChanged,
        wellsTemplateForm.onInputChanged,
        wellsOtherForm.onInputChanged
      ), 500).subscribe((_) => refreshResults());

    platesFormHost.appendChild(platesFormsResult.root);
    wellsFormHost.appendChild(wellsFormsResult.root);

    refreshResults();
  };
  const setTemplate = (template: PlateTemplate) => {
    plateTemplate = template;
    refreshUI();
  };

  const filterPanel = ui.divV([
    ui.divH([plateTemplateSelector.root, globalSearchInput.root]),
    ui.divH([
      ui.divV([ui.h2('Plates'), platesFormHost], {style: {flexGrow: '1'}}),
      ui.divV([ui.h2('Wells'), wellsFormHost], {style: {flexGrow: '1'}})
    ], {style: {width: '100%'}}),
  ], {classes: 'ui-panel'});

  filterPanel.style.overflowY = 'auto';
  filterPanel.style.padding = '0 12px';

  const splitter = ui.splitV(
    [filterPanel, resultsView.root],
    {style: {width: '100%', height: '100%'}},
    true
  );

  mainView.root.appendChild(splitter);

  setTemplate(plateTemplates[0]);
  return mainView;
}


export function searchPlatesView(): DG.View {
  return getSearchView(queryPlates, (grid) => {
    grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  });
}

export function searchWellsView(): DG.View {
  return getSearchView(queryWells, (grid) => {});
}
