/* eslint-disable max-len */
import {getPlateUniquePropertyValues, getWellUniquePropertyValues, initPlates, PlateQuery, PlateTemplate, PropertyCondition} from '../plates-crud';

import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NumericMatcher, StringInListMatcher} from '../numeric_matcher';
import * as rxjs from 'rxjs';
import {
  PlateProperty,
  queryPlates,
  queryWells,
  allProperties,
  plateTemplates
} from '../plates-crud';


type PropInput = DG.InputBase & { prop: PlateProperty };


function getSearchForm(
  allProperties: PlateProperty[],
  templateProperties: PlateProperty[],
  getUniqueValues: (prop: any) => string[]
): { root: HTMLElement, templateForm: DG.InputForm, otherForm: DG.InputForm } {
  const createInput = (prop: PlateProperty): PropInput | null => {
    if (prop.choices && prop.choices.length > 0) {
      let choicesList: string[];
      if (typeof prop.choices === 'string') {
        try {
          choicesList = JSON.parse(prop.choices as any);
        } catch {
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


function getSearchView(viewName: string,
  search: (query: PlateQuery) => Promise<DG.DataFrame>,
  onResults: (grid: DG.Grid) => void
): DG.View {
  const mainView = DG.View.create();
  const resultsView = DG.TableView.create(DG.DataFrame.create(0, viewName), false);
  resultsView.name = viewName;
  resultsView._onAdded();

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
      df.name = viewName;
      onResults(resultsView.grid);
    });
  };

  const refreshUI = () => {
    ui.empty(platesFormHost);
    ui.empty(wellsFormHost);

    const currentTemplatePlateProps = plateTemplate?.plateProperties?.map((p) => p as PlateProperty) ?? [];
    const currentTemplateWellProps = plateTemplate?.wellProperties?.map((p) => p as PlateProperty) ?? [];

    const allPlateProps = globalSearchInput.value ?
      Array.from(new Map(
        allProperties
          .filter((p) => p.scope === 'plate')
          .map((p) => [p.name, p])
      ).values()) :
      currentTemplatePlateProps;

    const allWellProps = globalSearchInput.value ?
      Array.from(new Map(
        allProperties
          .filter((p) => p.scope === 'well')
          .map((p) => [p.name, p])
      ).values()) :
      currentTemplateWellProps;

    const platesFormsResult = getSearchForm(allPlateProps, currentTemplatePlateProps, getPlateUniquePropertyValues);
    const wellsFormsResult = getSearchForm(allWellProps, currentTemplateWellProps, getWellUniquePropertyValues);

    platesTemplateForm = platesFormsResult.templateForm;
    platesOtherForm = platesFormsResult.otherForm;
    wellsTemplateForm = wellsFormsResult.templateForm;
    wellsOtherForm = wellsFormsResult.otherForm;

    [platesTemplateForm.root, platesOtherForm.root, wellsTemplateForm.root, wellsOtherForm.root].forEach((r) => r.style.width = '100%');

    DG.debounce(
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
    [filterPanel, resultsHost],
    {style: {width: '100%', height: '100%'}},
    true
  );

  mainView.root.appendChild(splitter);

  setTemplate(plateTemplates[0]);
  return mainView;
}


export function searchPlatesView(): DG.View {
  return getSearchView('Search Plates', queryPlates, (grid) => {
    grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await curves.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  });
}

export function searchWellsView(): DG.View {
  return getSearchView('Search Wells', queryWells, (grid) => {});
}
