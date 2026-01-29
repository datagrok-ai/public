/* eslint-disable max-len */
import {AnalysisCondition, AnalysisProperty, getPlateUniquePropertyValues, getWellUniquePropertyValues, initPlates, PlateQuery, PlateTemplate, PropertyCondition, getAnalysisRunGroups} from '../plates-crud';

import * as grok from 'datagrok-api/grok';
import {AnalysisManager} from '../../plate/analyses/analysis-manager';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {NumericMatcher, StringInListMatcher, StringMatcher} from '../matchers';
import * as rxjs from 'rxjs';
import {
  PlateProperty,
  queryPlates,
  queryWells,
  plateTemplates
} from '../plates-crud';


type PropInput = DG.InputBase & { prop: PlateProperty };
type AnalysisPropInput = DG.InputBase & { prop: AnalysisProperty };


function getSearchForm(
  allProperties: PlateProperty[],
  templateProperties: PlateProperty[],
  getUniqueValues: (prop: any) => string[]
): { root: HTMLElement, templateForm: DG.InputForm, otherForm: DG.InputForm } {
  const createInput = (prop: PlateProperty): PropInput | null => {
    let input: PropInput | null = null;

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
      input = ui.input.multiChoice(prop.name, {items: choicesList}) as PropInput;
    } else if (prop.type === DG.TYPE.STRING) {
      const uniqueValues = getUniqueValues(prop);
      if (uniqueValues.length > 0)
        input = ui.input.multiChoice(prop.name, {items: uniqueValues}) as PropInput;
      else
        input = ui.input.string(prop.name) as PropInput;
    } else if (prop.type === DG.TYPE.FLOAT || prop.type === DG.TYPE.INT) {
      input = ui.input.string(prop.name) as PropInput;
      input.addValidator((s) => NumericMatcher.parse(s) ? null : 'Invalid numerical criteria. Example: ">10"');
    } else if (prop.type === DG.TYPE.BOOL) {
      input = ui.input.bool(prop.name) as PropInput;
    }

    if (input)
      input.prop = prop;


    return input;
  };

  const templatePropNames = new Set(templateProperties.map((p) => p.name));

  const templateInputs = templateProperties.map(createInput).filter((i) => i) as PropInput[];
  const otherInputs = allProperties.filter((p) => !templatePropNames.has(p.name))
    .map(createInput).filter((i) => i) as PropInput[];

  const templateForm = DG.InputForm.forInputs(templateInputs);
  const otherForm = DG.InputForm.forInputs(otherInputs);

  const root = ui.divV([
    ...(templateInputs.length > 0 ? [
      ui.h3('Template Properties', {style: {marginTop: '0px', marginBottom: '2px'}}),
      templateForm.root
    ] : []),
    ...(otherInputs.length > 0 ? [otherForm.root] : [])
  ]);


  return {root, templateForm, otherForm};
}

function analysisFormToMatchers(form: DG.InputForm, analysisName: string): AnalysisCondition[] {
  const matchers: AnalysisCondition[] = [];
  if (!form || !analysisName) return matchers;

  for (const input of form.inputs) {
    if (input.value == null || input.value === '' || (Array.isArray(input.value) && input.value.length === 0)) continue;

    if (input.caption === 'Groups' && Array.isArray(input.value)) {
      matchers.push({
        analysisName: analysisName,
        group: input.value,
      });
      continue;
    }
    // Handle property conditions
    const propInput = input as AnalysisPropInput;
    if (propInput.prop && propInput.inputType === DG.InputType.Text && NumericMatcher.parse(propInput.value)) {
      matchers.push({
        property: propInput.prop,
        matcher: NumericMatcher.parse(propInput.value)!,
        analysisName: analysisName,
      });
    }
  }
  return matchers;
}

function searchFormToMatchers(form: DG.InputForm): PropertyCondition[] {
  const matchers: PropertyCondition[] = [];
  if (!form) return matchers;

  for (const input of form.inputs as PropInput[]) {
    const prop = input.prop;
    const value = input.value;

    if (value == null || value === '' || (Array.isArray(value) && value.length === 0))
      continue;

    if (input.inputType === DG.InputType.MultiChoice) {
      matchers.push({property: prop, matcher: new StringInListMatcher(value)});
    } else if (input.inputType === DG.InputType.Text) {
      if (prop.type === DG.TYPE.INT || prop.type === DG.TYPE.FLOAT) {
        const numericMatcher = NumericMatcher.parse(value);
        if (numericMatcher)
          matchers.push({property: prop, matcher: numericMatcher});
      } else if (prop.type === DG.TYPE.STRING) {
        matchers.push({property: prop, matcher: new StringMatcher(value)});
      }
    }
  }
  return matchers;
}


function getSearchView(viewName: string,
  search: (query: PlateQuery) => Promise<DG.DataFrame>,
  onResults: (grid: DG.Grid) => void
): DG.View {
  const mainView = DG.View.create();
  mainView.name = viewName;
  const resultsView = DG.TableView.create(DG.DataFrame.create(0, viewName), false);
  resultsView.name = viewName;
  resultsView._onAdded();

  let platesTemplateForm: DG.InputForm;
  let platesOtherForm: DG.InputForm;
  let wellsTemplateForm: DG.InputForm;
  let wellsOtherForm: DG.InputForm;
  let analysisForm: DG.InputForm;

  const platesFormHost = ui.div();
  const wellsFormHost = ui.div();
  const analysisFormHost = ui.div();

  let plateTemplate = plateTemplates[0];

  const plateTemplateSelector = ui.input.choice('Template', {
    items: plateTemplates.map((pt) => pt.name),
    value: plateTemplate?.name,
    onValueChanged: (v) => {
      const selectedTemplate = plateTemplates.find((pt) => pt.name === v);
      if (selectedTemplate)
        setTemplate(selectedTemplate);
    }
  });

  const analysisTypeSelector = ui.input.choice('Analysis', {
    items: ['None', ...AnalysisManager.instance.analyses.map((a) => a.friendlyName)],
    value: 'None',
    onValueChanged: () => refreshAnalysisUI(),
  });

  const refreshResults = () => {
    const selectedAnalysis = AnalysisManager.instance.byFriendlyName(analysisTypeSelector.value as string);
    const analysisName = selectedAnalysis ? selectedAnalysis.name : '';

    const analysisMatchers = analysisFormToMatchers(analysisForm, analysisName);

    // If an analysis is selected (analysisName is not 'None' or empty),
    // but no sub-filters are filled out (analysisMatchers is empty),
    // we must add a "base" condition to filter by the analysis type itself.
    if (analysisName && analysisMatchers.length === 0) {
      analysisMatchers.push({
        analysisName: analysisName
      } as AnalysisCondition);
    }

    const query: PlateQuery = {
      plateMatchers: [...searchFormToMatchers(platesTemplateForm), ...searchFormToMatchers(platesOtherForm)],
      wellMatchers: [...searchFormToMatchers(wellsTemplateForm), ...searchFormToMatchers(wellsOtherForm)],
      analysisMatchers: analysisMatchers,
    };

    search(query).then((df) => {
      resultsView.dataFrame = df;
      df.name = viewName;
      resultsView.grid.setOptions({allowEdit: false});
      onResults(resultsView.grid);
    });
  }; ;

  const refreshAnalysisUI = async () => {
    ui.empty(analysisFormHost);
    const selectedAnalysisName = analysisTypeSelector.value;

    if (selectedAnalysisName === 'None' || selectedAnalysisName === null) {
      analysisForm = DG.InputForm.forInputs([]);
      refreshResults();
      return;
    }

    const analysis = AnalysisManager.instance.byFriendlyName(selectedAnalysisName);
    if (!analysis) {
      console.error(`Analysis "${selectedAnalysisName}" not found.`);
      refreshResults();
      return;
    }

    const inputs: DG.InputBase[] = [];

    if (analysis.name === 'DRC' || analysis.name === 'Dose-Ratio') {
      const groups = await getAnalysisRunGroups(analysis.name);
      if (groups.length > 0) {
        const groupInput = ui.input.multiChoice('Groups', {
          items: groups,
          value: [],
        });
        inputs.push(groupInput);
      }
    }

    for (const prop of analysis.outputs) {
      if (prop.type === DG.TYPE.FLOAT || prop.type === DG.TYPE.INT) {
        const input = ui.input.string(prop.name) as AnalysisPropInput;
        input.prop = prop as AnalysisProperty;
        input.addValidator((s) => (s === '' || NumericMatcher.parse(s)) ? null : 'Invalid criteria. Example: ">10"');
        inputs.push(input);
      }
    }

    analysisForm = DG.InputForm.forInputs(inputs);
    analysisForm.root.style.width = '100%';
    analysisFormHost.appendChild(analysisForm.root);

    DG.debounce(analysisForm.onInputChanged, 500).subscribe((_) => refreshResults());

    refreshResults();
  };

  const refreshUI = () => {
    ui.empty(platesFormHost);
    ui.empty(wellsFormHost);

    const currentTemplatePlateProps = plateTemplate?.plateProperties?.map((p) => p as PlateProperty) ?? [];
    const currentTemplateWellProps = plateTemplate?.wellProperties?.map((p) => p as PlateProperty) ?? [];

    const allPlateProps = currentTemplatePlateProps;
    const allWellProps = currentTemplateWellProps;

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
    ui.divH([plateTemplateSelector.root]),
    ui.divH([
      ui.divV([ui.h2('Plates'), platesFormHost], {style: {flexGrow: '1T', marginRight: '10px'}}),
      ui.divV([ui.h2('Wells'), wellsFormHost], {style: {flexGrow: '1', marginLeft: '10px'}})
    ], {style: {width: '100%'}}),
    ui.divV([
      ui.h2('Analysis Results'),
      analysisTypeSelector.root,
      analysisFormHost
    ], {style: {marginTop: '10px'}})
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
  setTimeout(() => {
    refreshAnalysisUI();
  }, 1500);
  return mainView;
}


export function searchPlatesView(): DG.View {
  return getSearchView('Search Plates', queryPlates, (grid) => {
    grid.columns.add({gridColumnName: 'plate', cellType: 'Plate'})
      .onPrepareValueScript = `return (await plates.getPlateByBarcode(gridCell.tableRow.get('barcode'))).data;`;
  });
}

export function searchWellsView(): DG.View {
  return getSearchView('Search Wells', queryWells, (grid) => {});
}
