/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  DesirabilityMode, MpoScale, NumericalDesirability, PropertyDesirability, WEIGHTED_AGGREGATIONS_LIST,
  WeightedAggregation, getMpoScoreColumnName, isNumerical, setDesirabilityScale,
} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {MpoProfileCreateView} from '../mpo/mpo-create-profile';
import {MpoProfilesView} from '../mpo/mpo-profiles-view';
import {MpoProfileHandler} from '../mpo/mpo-profile-handler';
import {mpoProfileStore} from '../mpo/mpo-profile-store';
import {importProfile} from '../mpo/mpo-profile-actions';
import {computeMpo, MpoMethod, MpoProfileInfo} from '../mpo/utils';

const _funcs = new Map<string, DG.Func>();

function registerAi<T>(name: string, signature: string, description: string, resolve: (view: any) => T,
  run: (target: T, ...args: any[]) => Promise<any>): DG.Func {
  let f = _funcs.get(name);
  if (!f) {
    f = grok.functions.register({signature, isAsync: true, namespace: 'ChemMpoAI',
      tags: 'chem-mpo-ai-function', options: {description},
      run: async (view: any, ...args: any[]) => run(resolve(view), ...args)});
    _funcs.set(name, f);
  }
  return f;
}

function aiFunc(name: string, signature: string, description: string,
  run: (view: MpoProfileCreateView, ...args: any[]) => Promise<any>): DG.Func {
  return registerAi(name, signature, description, editorOf, run);
}

function profilesAiFunc(name: string, signature: string, description: string,
  run: (view: MpoProfilesView, ...args: any[]) => Promise<any>): DG.Func {
  return registerAi(name, signature, description, profilesViewOf, run);
}

const editorByView = new WeakMap<object, MpoProfileCreateView>();
const profilesByView = new WeakMap<object, MpoProfilesView>();

function registerView<T extends object>(map: WeakMap<object, T>, view: DG.ViewBase, target: T): DG.ViewBase {
  const wrapper = (view.dart ? DG.toJs(view.dart) : null) ?? view;
  map.set(view, target);
  map.set(wrapper, target);
  map.set(view.dart, target);
  return wrapper;
}

function resolveView<T>(map: WeakMap<object, T>, shellView: any): T | null {
  for (const candidate of [shellView, shellView?.dart, shellView?.jsView]) {
    const target = candidate != null ? map.get(candidate) : null;
    if (target)
      return target;
  }
  return null;
}

export function attachMpoProfileAi(editor: MpoProfileCreateView): void {
  const w = registerView(editorByView, editor.activeView, editor);
  (w as any).aiDescription = MPO_EDITOR_AI_DESCRIPTION;
  (w as any).getFunctions = () => [...mpoFunctions(), ...DG.View.prototype.getFunctions.call(w)];
}

export function attachMpoProfilesAi(profiles: MpoProfilesView): void {
  const w = registerView(profilesByView, profiles.view, profiles);
  (w as any).aiDescription = MPO_PROFILES_AI_DESCRIPTION;
  (w as any).getFunctions = () => mpoProfilesFunctions();
}

function editorOf(shellView: any): MpoProfileCreateView {
  const editor = resolveView(editorByView, shellView);
  if (!editor)
    throw new Error('The current view is not an MPO profile editor - open one from Browse | Chem | MPO profiles');
  return editor;
}

function profilesViewOf(shellView: any): MpoProfilesView {
  const view = resolveView(profilesByView, shellView);
  if (!view)
    throw new Error('The current view is not the MPO profile list - open it from Browse | Chem | MPO profiles');
  return view;
}

type FoundColumn = {col: DG.Column, error?: undefined} | {col?: undefined, error: string};
type FoundProfile = {profile: MpoProfileInfo, error?: undefined} | {profile?: undefined, error: string};

const NO_DATASET = 'The editor has no dataset - set one with setMpoDataset first';

function propertyError(view: MpoProfileCreateView, property: string): string | null {
  return view.profile.properties[property] ? null :
    `No property '${property}' - the profile scores ${Object.keys(view.profile.properties).join(', ')}`;
}

function columnOf(view: MpoProfileCreateView, column: string): FoundColumn {
  if (!view.df)
    return {error: NO_DATASET};
  const col = view.df.col(column);
  return col ? {col} : {error: `No column '${column}' in ${view.df.name} - call getMpoState for the ones it has`};
}

async function profileOf(name: string): Promise<FoundProfile> {
  const profiles = await mpoProfileStore.ensureLoaded();
  const wanted = name?.trim().toLowerCase();
  const profile = profiles.find((p) => p.name.toLowerCase() === wanted);
  return profile ? {profile} : {error: `No profile '${name}' - call listMpoProfiles for the ones there are`};
}

function aggregationOf(view: MpoProfileCreateView): WeightedAggregation {
  return view.editor.aggregationInput.value ?? view.profile.aggregation ?? 'Average';
}

function scoreColumn(view: MpoProfileCreateView): DG.Column | null {
  return view.df?.col(getMpoScoreColumnName(view.profile.name)) ?? null;
}

function reportedProperty(prop: PropertyDesirability): object {
  const reported: any = {...prop};
  delete reported.freeformLine;
  if (isNumerical(prop) && prop.mode !== DesirabilityMode.Freeform)
    delete reported.line;
  return reported;
}

function state(view: MpoProfileCreateView): object {
  const df = view.df;
  const score = scoreColumn(view);
  return {
    profile: {...view.profile, properties: Object.fromEntries(
      Object.entries(view.profile.properties).map(([name, prop]) => [name, reportedProperty(prop)]))},
    columnMapping: view.editor.columnMapping,
    aggregation: aggregationOf(view),
    method: view.showMethod ? (view.isManualMode ? MpoMethod.Manual : MpoMethod.DataDriven) : null,
    modified: view.isModified,
    dataset: df ? {name: df.name, rowCount: df.rowCount, columns: df.columns.names()} : null,
    unmapped: Object.keys(view.profile.properties)
      .filter((p) => !df?.col(view.editor.columnMapping[p] ?? p)),
    score: score ? {column: score.name, min: score.stats.min, max: score.stats.max, avg: score.stats.avg} : null,
  };
}

function changed(view: MpoProfileCreateView): object {
  return {success: true, ...state(view)};
}

const MPO_EDITOR_AI_DESCRIPTION = 'MPO profile editor - a multi-parameter optimization profile: the scored properties, ' +
  'each with a weight and a desirability curve, previewed against a dataset. Act on it through the view ' +
  'functions (search list_view_functions with "mpo"): getMpoState (call first - it gives the property names, ' +
  'the column each is mapped to, and the desirability fields setMpoDesirability takes back), then ' +
  'addMpoProperty, removeMpoProperty, renameMpoProperty, setMpoPropertyColumn, setMpoPropertyWeight and ' +
  'setMpoDesirability for the properties; setMpoProfileName, setMpoProfileDescription, setMpoDataset, ' +
  'setMpoAggregation and setMpoMethod for the profile; saveMpoProfile, resetMpoProfile and applyMpoProfile ' +
  'for the rest.';

export function mpoFunctions(): DG.Func[] {
  return [
    aiFunc('getMpoState', 'dynamic getMpoState(view view)',
      'The profile as it stands: every property with its weight, desirability function and missing-value handling - including the full point list of a freeform curve - plus the column each is mapped to, the aggregation, the method, the dataset columns, and the score column stats. Call this first: the other functions address a property by its name here, and setMpoDesirability takes the same desirability fields back',
      async (view) => state(view)),

    aiFunc('addMpoProperty', 'dynamic addMpoProperty(view view, string property, string column)',
      'Start scoring a new property, named by property. column is optional - it maps the property to a dataset column, which seeds the curve from that column: its range for a numerical column, its categories for a categorical one. The property is added either way, so pass column only when a dataset is loaded and you know which column holds the values. Shape the curve afterwards with setMpoDesirability',
      async (view, property: string, column?: string) => {
        const name = property?.trim();
        if (!name)
          return {success: false, error: 'Pass property - the name of the property to add'};
        if (view.profile.properties[name])
          return {success: false, error: `'${name}' is already scored - change it with setMpoDesirability, or pick another name`};

        view.editor.addProperty(name);
        const found = column == null ? null : columnOf(view, column);
        if (found?.col)
          view.editor.setPropertyColumn(name, found.col);
        return found?.error ? {...changed(view), note: `'${name}' was added but left unmapped: ${found.error}`} : changed(view);
      }),

    aiFunc('removeMpoProperty', 'dynamic removeMpoProperty(view view, string property)',
      'Stop scoring a property. The profile keeps the rest, and the score is recomputed without it',
      async (view, property: string) => {
        const error = propertyError(view, property);
        if (error)
          return {success: false, error};
        view.editor.deleteProperty(property);
        return changed(view);
      }),

    aiFunc('renameMpoProperty', 'dynamic renameMpoProperty(view view, string property, string newName)',
      'Rename a scored property, keeping its curve, weight and column mapping',
      async (view, property: string, newName: string) => {
        const error = propertyError(view, property);
        if (error)
          return {success: false, error};
        const name = newName?.trim();
        if (!name)
          return {success: false, error: 'Pass newName - the name to give the property'};
        if (!view.editor.renameProperty(property, name))
          return {success: false, error: `Could not rename '${property}' to '${name}' - a property with that name already exists`};
        return changed(view);
      }),

    aiFunc('setMpoPropertyColumn', 'dynamic setMpoPropertyColumn(view view, string property, string column)',
      'Map a property to the dataset column holding its values; the curve switches to match the column type, and a numerical range the user has not pinned is refitted to the column. Leave column out to unmap the property, which stops it scoring until it is mapped again',
      async (view, property: string, column?: string) => {
        const error = propertyError(view, property);
        if (error)
          return {success: false, error};

        const found = column == null ? null : columnOf(view, column);
        if (found?.error)
          return {success: false, error: found.error};

        view.editor.setPropertyColumn(property, found?.col ?? null);
        return changed(view);
      }),

    aiFunc('setMpoPropertyWeight', 'dynamic setMpoPropertyWeight(view view, string property, double weight)',
      'How much a property counts towards the score, 0 to 1, relative to the other weights',
      async (view, property: string, weight: number) => {
        const error = propertyError(view, property);
        if (error)
          return {success: false, error};
        view.editor.updateProperty(property, {weight: Math.max(0, Math.min(1, weight))});
        return changed(view);
      }),

    aiFunc('setMpoDesirability', 'dynamic setMpoDesirability(view view, string property, dynamic desirability)',
      'Shape the curve mapping a property\'s values to a 0-1 desirability. desirability patches the property\'s own fields, exactly as getMpoState reports them: mode (freeform, gaussian or sigmoid) with mean and sigma for a gaussian peak or x0 and k for a sigmoid rise, min and max for the range, line as the complete list of [value, desirability] pairs for a freeform curve - getMpoState reports the points it has now, so send those back with your change applied rather than only the new ones - inverted true to flip the curve (lower is better), scale log for values spanning orders of magnitude - the curve keeps its shape across the flip - categories as {name, desirability} for a categorical property (cover every category, uncovered values score the whole row null), and missingValues {strategy} of exclude, skip, or default with a score',
      async (view, property: string, desirability: Partial<PropertyDesirability>) => {
        const error = propertyError(view, property);
        if (error)
          return {success: false, error};
        if (!desirability)
          return {success: false, error: 'Pass desirability - the fields of the curve to change, as getMpoState reports them'};

        const prop = view.profile.properties[property];
        const patch = {...desirability} as Partial<NumericalDesirability>;
        if (isNumerical(prop) && patch.scale != null) {
          setDesirabilityScale(prop, patch.scale === MpoScale.Log);
          delete patch.scale;
        }
        if (isNumerical(prop) && (patch.min != null || patch.max != null))
          patch.rangeUserSet = true;

        view.editor.updateProperty(property, patch);
        return changed(view);
      }),

    aiFunc('setMpoProfileName', 'dynamic setMpoProfileName(view view, string name)',
      'Rename the profile. The score column is renamed along with it; saveMpoProfile still writes to the file the profile came from',
      async (view, name: string) => {
        const newName = name?.trim();
        if (!newName)
          return {success: false, error: 'Pass name - what to call the profile'};
        view.setProfileName(newName);
        return changed(view);
      }),

    aiFunc('setMpoProfileDescription', 'dynamic setMpoProfileDescription(view view, string description)',
      'Describe what the profile optimizes for - shown in the profile list',
      async (view, description: string) => {
        view.setProfileDescription(description ?? '');
        return changed(view);
      }),

    aiFunc('setMpoDataset', 'dynamic setMpoDataset(view view, dataframe dataset)',
      'The table the curves are previewed and scored against. dataset is a table open in the workspace, not a file - open the file first if it is not open yet. Switching it re-maps the properties onto the new columns, and asks the user what to do with unsaved curve changes',
      async (view, dataset: DG.DataFrame) => {
        view.setDataset(dataset ?? null);
        return changed(view);
      }),

    aiFunc('setMpoAggregation', 'dynamic setMpoAggregation(view view, string aggregation)',
      `How the per-property desirabilities combine into one score. One of ${WEIGHTED_AGGREGATIONS_LIST.join(', ')}: Average and Geomean normalize by the total weight, while Sum is unbounded and Product, Min and Max apply weights as d^w, so changing it changes what the weights mean`,
      async (view, aggregation: string) => {
        const wanted = aggregation?.toLowerCase();
        const agg = WEIGHTED_AGGREGATIONS_LIST.find((a) => a.toLowerCase() === wanted);
        if (!agg)
          return {success: false, error: `aggregation must be one of: ${WEIGHTED_AGGREGATIONS_LIST.join(', ')}`};
        view.editor.aggregationInput.value = agg;
        return changed(view);
      }),

    aiFunc('setMpoMethod', 'dynamic setMpoMethod(view view, string method)',
      `How the profile is built: ${MpoMethod.Manual} to shape the curves by hand, or ${MpoMethod.DataDriven} to train them from the dataset - the properties that best separate the compounds labelled preferred from the rest, each with a fitted curve and a weight, plus a ROC curve and confusion matrix. Training needs the EDA package and a labelled column, takes a while, and sets the current curves aside, so only switch to it when the user asks for a data-driven profile`,
      async (view, method: string) => {
        const wanted = method?.toLowerCase();
        const value = [MpoMethod.Manual, MpoMethod.DataDriven].find((m) => m.toLowerCase() === wanted);
        if (!value)
          return {success: false, error: `method must be one of: ${MpoMethod.Manual}, ${MpoMethod.DataDriven}`};
        if (!view.showMethod)
          return {success: false, error: 'The method is fixed while editing a saved profile - it only applies while creating one'};
        if (value === MpoMethod.DataDriven && !view.df)
          return {success: false, error: NO_DATASET};

        await view.setMethod(value);
        return changed(view);
      }),

    aiFunc('saveMpoProfile', 'dynamic saveMpoProfile(view view, string name)',
      'Save the profile. Pass name to save it under a different name, which leaves the original untouched - that is how a profile is cloned. Without a name it overwrites the profile being edited',
      async (view, name?: string) => {
        const newName = name?.trim();
        if (newName && newName !== view.profile.name) {
          view.setProfileName(newName);
          view.saved = null;
        }
        if (!await view.save())
          return {success: false, error: `Profile '${view.profile.name}' was not saved`};
        return {success: true, name: view.profile.name, profileId: view.saved?.id};
      }),

    aiFunc('resetMpoProfile', 'dynamic resetMpoProfile(view view)',
      'Discard every change made since the profile was last saved',
      async (view) => {
        await view.resetProfile();
        return changed(view);
      }),

    aiFunc('applyMpoProfile', 'dynamic applyMpoProfile(view view, bool desirabilityColumns)',
      'Score the dataset with the profile as it stands and add the score column to the table as a recorded transform - the editor only previews it otherwise. desirabilityColumns also adds one 0-1 column per property, showing what each contributes to a row\'s score',
      async (view, desirabilityColumns?: boolean) => {
        if (!view.df)
          return {success: false, error: NO_DATASET};
        const added = await computeMpo(view.df, view.profile, view.editor.columnMapping, aggregationOf(view),
          false, false, !!desirabilityColumns);
        if (added.length === 0)
          return {success: false, error: 'Nothing was scored - no profile property resolves to a column of this table'};
        return {success: true, scoreColumn: added[0], ...state(view)};
      }),
  ];
}

const MPO_PROFILES_AI_DESCRIPTION = 'MPO profile list - the saved multi-parameter optimization profiles, each ' +
  'scoring a set of properties by desirability. Act on it through the view functions (search ' +
  'list_view_functions with "mpo"): listMpoProfiles (call first), previewMpoProfile to show one in the context ' +
  'panel, openMpoProfile / createMpoProfile / cloneMpoProfile to open the editor, which reports its own set of ' +
  'functions once open, deleteMpoProfile, and importMpoProfile to load one from a JSON file on the server.';

export function mpoProfilesFunctions(): DG.Func[] {
  return [
    profilesAiFunc('listMpoProfiles', 'dynamic listMpoProfiles(view view)',
      'Every saved MPO profile - its name, description, the properties it scores and how they aggregate. Call this first; the other functions address a profile by the name reported here',
      async () => {
        const profiles = await mpoProfileStore.ensureLoaded();
        return profiles.map((p) => ({
          name: p.name,
          description: p.description,
          id: p.id,
          aggregation: p.aggregation ?? null,
          properties: Object.keys(p.properties ?? {}),
        }));
      }),

    profilesAiFunc('previewMpoProfile', 'dynamic previewMpoProfile(view view, string profile)',
      'Show a profile in the context panel - its properties with their curves and weights - without opening it for editing. profile is the name listMpoProfiles reports',
      async (view, profile: string) => {
        const {profile: found, error} = await profileOf(profile);
        if (!found)
          return {success: false, error};
        view.preview(found);
        return {success: true, name: found.name, properties: Object.keys(found.properties ?? {})};
      }),

    profilesAiFunc('openMpoProfile', 'dynamic openMpoProfile(view view, string profile)',
      'Open a saved profile in the editor, where its properties, curves, weights and column mapping can be changed. profile is the name listMpoProfiles reports. The editor reports its own functions - call list_view_functions again once it is open',
      async (_view, profile: string) => {
        const {profile: found, error} = await profileOf(profile);
        if (!found)
          return {success: false, error};
        MpoProfileHandler.edit(found);
        return {success: true, name: found.name, note: 'The profile editor is open'};
      }),

    profilesAiFunc('createMpoProfile', 'dynamic createMpoProfile(view view, string name)',
      'Open the editor on a new profile, called name if one is given. It starts with no properties - add them there with addMpoProperty - and nothing is saved until saveMpoProfile is called',
      async (view, name?: string) => {
        const created = view.openCreateProfile(name?.trim());
        return {success: true, name: created.profile.name, note: 'The profile editor is open on the new profile'};
      }),

    profilesAiFunc('cloneMpoProfile', 'dynamic cloneMpoProfile(view view, string profile)',
      'Open a copy of a saved profile in the editor, leaving the original untouched. profile is the name listMpoProfiles reports. The copy is not saved until saveMpoProfile is called there',
      async (_view, profile: string) => {
        const {profile: found, error} = await profileOf(profile);
        if (!found)
          return {success: false, error};
        MpoProfileHandler.clone(found);
        return {success: true, clonedFrom: found.name, note: 'The profile editor is open on the copy'};
      }),

    profilesAiFunc('deleteMpoProfile', 'dynamic deleteMpoProfile(view view, string profile)',
      'Delete a saved profile, named by profile. The user is asked to confirm, and nothing is deleted until they do',
      async (_view, profile: string) => {
        const {profile: found, error} = await profileOf(profile);
        if (!found)
          return {success: false, error};
        MpoProfileHandler.delete(found);
        return {success: true, name: found.name, note: 'Asked the user to confirm the deletion'};
      }),

    profilesAiFunc('importMpoProfile', 'dynamic importMpoProfile(view view, string file)',
      'Add a profile from a JSON file on the server, given its path, e.g. System:AppData/Home/lipinski.json. A name already taken asks the user whether to replace that profile or keep both',
      async (_view, file: string) => {
        const path = file?.trim();
        if (!path)
          return {success: false, error: 'Pass file - the path of the .json profile to import'};
        const name = path.split(/[\\/]/).pop()!.replace(/\.json$/i, '');
        if (!await importProfile(await grok.dapi.files.readAsText(path), name))
          return {success: false, error: `'${path}' was not imported - it is not a valid MPO profile, or the user cancelled`};
        return {success: true, imported: path};
      }),
  ];
}
