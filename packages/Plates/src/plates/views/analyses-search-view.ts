/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {AnalysisManager} from '../../plate/analyses/analysis-manager';
import {AnalysisProperty, AnalysisQuery, getAnalysisRunGroups, queryAnalysesGeneric, PropertyCondition, allProperties, queryAnalyses} from '../plates-crud';
import {NumericMatcher} from '../matchers';
import './analyses-search-view.css';

type AnalysisPropInput = DG.InputBase & { prop: AnalysisProperty };

function analysisFormToQuery(form: DG.InputForm, analysisName: string): Omit<AnalysisQuery, 'analysisName'> {
  const propertyMatchers: PropertyCondition[] = [];
  let group: string | string[] | undefined;

  if (!form || !analysisName) return {propertyMatchers, group};

  for (const input of form.inputs) {
    if (input.value == null || input.value === '' || (Array.isArray(input.value) && input.value.length === 0))
      continue;

    if (input.caption === 'Groups' && Array.isArray(input.value) && input.value.length > 0) {
      group = input.value;
      continue;
    }
    const propInput = input as AnalysisPropInput;
    const analysisProp = propInput.prop;

    if (analysisProp && (analysisProp.type === DG.TYPE.FLOAT || analysisProp.type === DG.TYPE.INT)) {
      const matcher = NumericMatcher.parse(propInput.value);
      if (matcher) {
        const fullPlateProperty = allProperties.find((p) => p.name === analysisProp.name);

        if (fullPlateProperty) {
          propertyMatchers.push({
            property: fullPlateProperty,
            matcher: matcher,
          });
        } else {
          grok.shell.warning(`Could not find property details for "${analysisProp.name}". This filter will be ignored.`);
        }
      }
    }
  }
  return {propertyMatchers, group};
}

function getAnalysesSearchView(
  viewName: string,
  search: (query: AnalysisQuery) => Promise<DG.DataFrame>,
  onResults: (grid: DG.Grid, analysisName: string) => void
): DG.View {
  const mainView = DG.View.create();
  mainView.name = viewName;
  mainView.root.classList.add('assay-analyses-search-view');
  const initialDf = DG.DataFrame.create(0);
  const resultsView = DG.TableView.create(initialDf);
  resultsView.root.classList.add('assay-analyses-search-view__results-container');
  let analysisForm: DG.InputForm;

  const analysisFormHost = ui.div([], 'assay-analyses-search-view__form-host');

  const refreshResults = () => {
    const selectedAnalysis = AnalysisManager.instance.byFriendlyName(analysisTypeSelector.value as string);
    if (!selectedAnalysis) {
      resultsView.dataFrame = DG.DataFrame.create(0);
      return;
    }

    const queryBase = analysisFormToQuery(analysisForm, selectedAnalysis.name);
    const query: AnalysisQuery = {
      ...queryBase,
      analysisName: selectedAnalysis.name,
    };

    const searchPromise = selectedAnalysis.queryResults ?
      selectedAnalysis.queryResults(query) :
      queryAnalyses(query);

    searchPromise.then((df) => {
      if (df) {
        resultsView.dataFrame = df;
        df.name = viewName;
        onResults(resultsView.grid, selectedAnalysis.name);
        resultsView.grid.invalidate();
        setTimeout(() => {
          resultsView.grid.invalidate();
        }, 100);
      }
    }).catch((error) => {
      console.error('Error loading results:', error);
      grok.shell.error(`Failed to load results: ${error.message}`);
    });
  };

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

        groupInput.root.classList.add('assay-analyses-search-view__group-input');
        inputs.push(groupInput);
      }
    }
    const searchableProps = analysis.getSearchableProperties ?
      analysis.getSearchableProperties() :
      analysis.outputs;
    for (const prop of searchableProps) {
      if (prop.type === DG.TYPE.FLOAT || prop.type === DG.TYPE.INT) {
        const input = ui.input.string(prop.name) as AnalysisPropInput;
        input.prop = prop as AnalysisProperty;
        input.addValidator((s) =>
          (s === '' || NumericMatcher.parse(s)) ? null : 'Invalid criteria. Example: ">10"'
        );
        inputs.push(input);
      }
    }

    analysisForm = DG.InputForm.forInputs(inputs);
    analysisForm.root.style.width = '100%';
    analysisFormHost.appendChild(analysisForm.root);

    DG.debounce(analysisForm.onInputChanged, 500).subscribe((_) => refreshResults());

    refreshResults();
  };

  const analysisTypeSelector = ui.input.choice('Analysis', {
    items: ['None', ...AnalysisManager.instance.analyses.map((a) => a.friendlyName)],
    value: 'None',
    onValueChanged: refreshAnalysisUI,
  });
  analysisTypeSelector.root.classList.add('assay-analyses-search-view__analysis-type-selector');

  const filterPanel = ui.divV([
    ui.h2('Analysis Filters'),
    analysisTypeSelector.root,
    analysisFormHost,
  ], 'assay-analyses-search-view__filter-panel');
  const mainContainer = ui.divV([
    filterPanel,
    resultsView.root
  ], 'assay-analyses-search-view__main-container');

  mainView.root.appendChild(mainContainer);
  refreshAnalysisUI();
  return mainView;
}
export function searchAnalysesView(): DG.View {
  return getAnalysesSearchView('Search Analyses', queryAnalysesGeneric, (grid, analysisName) => {
    grid.setOptions({allowEdit: true});

    const plateIdCol = grid.col('plate_id');
    if (plateIdCol) plateIdCol.visible = false;

    const runIdCol = grid.col('run_id');
    if (runIdCol) runIdCol.visible = false;

    const groupCol = grid.col('group_combination');
    if (groupCol) groupCol.visible = false;

    if (grid.dataFrame.columns.byName('Curve')) {
      const curveCol = grid.dataFrame.getCol('Curve');
      curveCol.semType = 'fit';
      curveCol.setTag(DG.TAGS.CELL_RENDERER, 'fit');

      const s = DG.debounce(grid.onAfterDrawContent, 300).subscribe(() => {
        s.unsubscribe();
        const visibleCurveCol = grid.col('Curve');
        if (visibleCurveCol) {
          visibleCurveCol.width = grid.root.clientWidth / 3;
          grid.props.rowHeight = 250;


          grid.invalidate();
        }
      });
    }
  });
}
