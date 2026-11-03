import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {DEFAULT_FLOAT_FORMAT} from '@datagrok-libraries/webcomponents-vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {tryParseJson} from '@datagrok-libraries/utils/src/string-utils';
import {useObservable} from '@vueuse/rxjs';
import {startWith, map} from 'rxjs/operators';
import type {ScalarsSection, ScalarState} from './ScalarsPanel';

export interface ScalarsState {
  type: 'scalars',
  scalarsData: ScalarState[],
  sections?: ScalarsSection[],
}

type OutputCategoryGroupSpec = ReadonlyArray<string | Record<string, OutputCategoryGroupSpec>>;
type OutputCategoryGroups = Record<string, OutputCategoryGroupSpec>;

function parseOutputCategoryGroups(func: DG.Func): OutputCategoryGroups | undefined {
  const raw = func?.options?.['outputCategoryGroups'];
  if (!raw) return undefined;
  const parsed = typeof raw === 'string' ? tryParseJson(raw, {warning: 'outputCategoryGroups'}) : raw;
  return (parsed && typeof parsed === 'object' && !Array.isArray(parsed)) ? parsed : undefined;
}

function applyOutputCategoryGroups(outputs: TabContent, spec: OutputCategoryGroups) {
  const consume = (items: OutputCategoryGroupSpec, depth: number, sink: ScalarsSection[]) => {
    for (const item of items) {
      if (typeof item === 'string') {
        const entry = outputs.get(item);
        if (!entry || entry.type !== 'scalars' || entry.scalarsData.length === 0) continue;
        sink.push({label: item, indent: depth, scalarsData: entry.scalarsData});
        outputs.delete(item);
      } else if (item && typeof item === 'object' && !Array.isArray(item)) {
        for (const [label, nested] of Object.entries(item)) {
          if (!Array.isArray(nested)) continue;
          const start = sink.length;
          consume(nested, depth + 1, sink);
          if (sink.length > start)
            sink.splice(start, 0, {label, indent: depth, scalarsData: []});
        }
      }
    }
  };

  for (const [superLabel, items] of Object.entries(spec)) {
    if (!Array.isArray(items)) continue;
    const sections: ScalarsSection[] = [];
    consume(items, 0, sections);
    if (sections.length === 0) continue;
    const flat = sections.flatMap((s) => s.scalarsData);
    outputs.set(superLabel, {type: 'scalars', scalarsData: flat, sections});
  }
}

export interface DataFrameState {
  name: string,
  df: Vue.ShallowRef<DG.DataFrame>,
  type: 'dataframe',
  config: Record<string, any>,
}

export interface GraphicsState {
  name: string,
  image: Vue.ShallowRef<string | null>,
  type: 'graphics',
}

export type TabContent = Map<string, ScalarsState | DataFrameState | GraphicsState>;

export interface RenderStateItem {
  tabLabel: string;
  tabContent: ScalarsState | DataFrameState | GraphicsState;
  isInput: boolean;
}

export const getEmptyTabToProperties = () => ({
  inputs: new Map() as TabContent,
  outputs: new Map() as TabContent,
});

const isEmptyScalar = (v: any) =>
  v == null || v === '' || v === DG.FLOAT_NULL || v === DG.INT_NULL;

const isEmptyDataFrame = (df: any) =>
  !df || (df instanceof DG.DataFrame && df.rowCount === 0);

const isEmptyGraphics = (v: any) => v == null || v === '';

const isHiddenOutput = (prop: DG.Property) => !!Utils.parseJsonOption(prop.options, 'hidden');

const getScalarContent = (funcCall: DG.FuncCall, prop: DG.Property) => {
  if (isHiddenOutput(prop))
    return;
  const scalarValue = funcCall.outputs[prop.name];
  let formattedScalarValue = scalarValue;

  if (prop.propertyType === DG.TYPE.FLOAT && scalarValue != null) {
    if (prop.options.format)
      formattedScalarValue = DG.format(scalarValue, prop.options.format);
    else if (prop.options.precision)
      formattedScalarValue = scalarValue.toPrecision(prop.options.precision);
    else
      formattedScalarValue = DG.format(scalarValue, DEFAULT_FLOAT_FORMAT);
  } else if (typeof scalarValue === 'boolean')
    formattedScalarValue = String(scalarValue);
  const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

  return [scalarValue, formattedScalarValue, units] as const;
};

// Stable identity per content so Vue skips re-applying `options` (and a redundant Dart setOptions) on rebuilds.
const viewerConfigIdentity = new Map<string, Record<string, string | boolean>>();
const stabilizeViewerConfig = (config: Record<string, string | boolean>) => {
  const key = JSON.stringify(config);
  const cached = viewerConfigIdentity.get(key);
  if (cached)
    return cached;
  viewerConfigIdentity.set(key, config);
  return config;
};

export const tabToProperties = (fc: DG.FuncCall) => {
  const tabsToProps = getEmptyTabToProperties();
  const hideEmpty = !Utils.getFeature(Utils.getFeatures(fc.func), 'show-empty-outputs', false);

  // Dataframes (one panel per viewer annotation) and graphics (a single image panel)
  const processVisual = (prop: DG.Property, isOutput: boolean) => {
    const isGraphics = prop.propertyType === DG.TYPE.GRAPHICS;
    let viewers = isGraphics ? [{}] : Utils.getPropViewers(prop).config;
    // Outputs without a viewer annotation (e.g. queries) get a plain grid by default
    const isDefaultGrid = !isGraphics && viewers.length === 0 && isOutput;
    if (isDefaultGrid)
      viewers = [{type: DG.VIEWER.GRID}];
    if (viewers.length === 0) return;
    const isEmpty = isGraphics ? isEmptyGraphics : isEmptyDataFrame;
    if (hideEmpty && isOutput && isEmpty(fc.outputs[prop.name])) return;

    viewers.forEach((viewer) => {
      const blockTitle = viewer.title ?? prop.options['caption'] ?? prop.name ?? ' ';
      const nameWithViewer = (isGraphics || isDefaultGrid) ? `${blockTitle}` : `${blockTitle} / ${viewer['type']}`;

      const tabLabel = prop.category === 'Misc' ?
        nameWithViewer: `${prop.category}: ${nameWithViewer}`;

      const name = prop.name;
      const source = isOutput ? fc.outputParams : fc.inputParams;
      const changes$ = source[name].onChanged.pipe(
        startWith(null),
        map(() => source[name].value ? Vue.markRaw(source[name].value) : null),
      );
      const value = useObservable(changes$);
      const content: DataFrameState | GraphicsState = isGraphics ?
        {type: 'graphics', name, image: value} :
        {type: 'dataframe', name, df: value, config: stabilizeViewerConfig(viewer)};
      if (isOutput)
        tabsToProps.outputs.set(tabLabel, content);
      else
        tabsToProps.inputs.set(tabLabel, content);
    });
    return;
  };

  [...fc.inputParams.values()].forEach(({ property }) => {
    if (property.propertyType === DG.TYPE.DATA_FRAME) processVisual(property, false);
  });

  [...fc.outputParams.values()].forEach(({property}) => {
    if (property.propertyType === DG.TYPE.DATA_FRAME || property.propertyType === DG.TYPE.GRAPHICS) {
      if (!isHiddenOutput(property))
        processVisual(property, true);
      return;
    }
    const category = property.category === 'Misc' ? 'Output': property.category;

    const categoryProps = tabsToProps.outputs.get(category);
    const content = getScalarContent(fc, property);
    if (!content)
      return;
    const [rawValue, formattedValue, units] = content;
    if (hideEmpty && isEmptyScalar(rawValue))
      return;
    const scalarProp = {name: property.name, friendlyName: property.caption || property.name, rawValue, formattedValue, units};
    if (categoryProps && categoryProps.type === 'scalars')
      categoryProps.scalarsData.push(scalarProp);
    else
      tabsToProps.outputs.set(category, {type: 'scalars', scalarsData: [scalarProp]});
  });

  const groupsSpec = parseOutputCategoryGroups(fc.func);
  if (groupsSpec)
    applyOutputCategoryGroups(tabsToProps.outputs, groupsSpec);

  return tabsToProps;
};
