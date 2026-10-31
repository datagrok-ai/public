/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import {
  CUSTOM_VIEW_OPTION, deepCopy, getRunAuthorId, getRunCreated, getRunName, isIncomplete, CREATED_BY_OPTION,
  CREATED_AT_OPTION, NESTED_ENTITY_OPTION,
} from '../../shared-utils/utils';
import dayjs from 'dayjs';

/** Name, description and tags of a saved run. */
export type RunMetadata = {
  name: string,
  description?: string,
  tags?: string[],
};

/** Where a run is saved: a space, or My runs when spaceId is not set. Runs in `replace` lose their space link. */
export type SaveTarget = {
  spaceId?: string,
  replace?: DG.FuncCall[],
};

export type NameConflict = {
  runs: DG.FuncCall[],
  other?: DG.Entity,
};

/** Fired on the platform's event bus after a run is saved on its own, with the spaces whose content changed. */
export const RUN_SAVED_EVENT = 'compute-run-saved';

export const NOT_DELETED = `((options.isDeleted != "true") or (options.isDeleted = null))`;

const normalizeName = (name: string) => name.trim().toLowerCase();

export namespace historyUtils {
  const groupsCache = new DG.LruCache();

  export function setRunMetadata(call: DG.FuncCall, metadata: RunMetadata) {
    call.friendlyName = metadata.name.trim();
    // the legacy title is read first, by Compute (v1) too, so it carries the same name
    call.options['title'] = call.friendlyName;
    // what isn't given stays as it is
    const current = getRunMetadata(call);
    call.options['description'] = metadata.description ?? current.description;
    call.options['tags'] = metadata.tags ?? current.tags;
  }

  export function getRunMetadata(call: DG.FuncCall): Required<RunMetadata> {
    const tags = call.options?.['tags'];
    return {
      name: getRunName(call),
      description: call.options?.['description'] ?? '',
      tags: Array.isArray(tags) ? tags : [],
    };
  }

  // a metadata save makes the saving user the platform author, so a run saved before the author option
  // existed records its original author and time first
  async function keepRunAuthor(call: DG.FuncCall) {
    if (call.options[CREATED_BY_OPTION])
      return;
    const source = getRunAuthorId(call) ? call : await shallowLoadRun(call.id);
    call.options[CREATED_BY_OPTION] = getRunAuthorId(source);
    call.options[CREATED_AT_OPTION] = getRunCreated(source)?.toISOString();
  }

  /** Edit details is allowed for the run's author, or for a user with Edit on the space the run is shown in. */
  export function canEditRunDetails(authorId: string | undefined, canEditSpace: boolean): boolean {
    return (!!authorId && authorId === DG.User.current().id) || canEditSpace;
  }

  /**
   * Shallow loads a FuncCall with a specified ID.
   * @param callToSave FuncCall to save
   */
  export async function shallowLoadRun(funcCallId: string) {
    const pulledRun = await grok.dapi.functions.calls.allPackageVersions()
      .include('session.user,func.package, inputs, outputs').find(funcCallId);
    return pulledRun;
  }

  /**
   * Loads a FuncCall with a specified ID. By default, also loads its' inputs/outputs and author.
   * FuncCall is loaded with internal TableInfo structs instead of DG.Dataframe-s.
   * Thus, we should load them separately, and it is time-consuming.
   * @returns Requested FuncCall
   */
  export async function loadRun(funcCallId: string) {
    const pulledRun = await grok.dapi.functions.calls.allPackageVersions()
      .include('session.user,func.package, inputs, outputs').find(funcCallId);

    const fileInputs = wu(pulledRun.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.FILE &&
        !!pulledRun.inputs[input.name]);

    await Promise.all(fileInputs
      .map(async (input) => {
        const {id, name} = JSON.parse(pulledRun.inputs[input.name]);
        const bytes = await grok.dapi.files.readAsBytes(id);
        const fileInfo = DG.FileInfo.fromBytes(name, bytes);
        fileInfo.id = id;
        pulledRun.inputs[input.name] = fileInfo;

        return Promise.resolve();
      }),
    );

    const dfOutputs = wu(pulledRun.outputParams.values() as DG.FuncCallParam[])
      .filter((output) =>
        output.property.propertyType === DG.TYPE.DATA_FRAME &&
        !!pulledRun.outputs[output.name],
      );
    await Promise.all(dfOutputs.map(async (output) => {
      pulledRun.outputs[output.name] = await grok.dapi.tables.getTable(pulledRun.outputs[output.name]);
      return Promise.resolve();
    }));

    const dfInputs = wu(pulledRun.inputParams.values() as DG.FuncCallParam[])
      .filter((input) =>
        input.property.propertyType === DG.TYPE.DATA_FRAME &&
        !!pulledRun.inputs[input.name],
      );
    await Promise.all(dfInputs.map(async (input) => {
      pulledRun.inputs[input.name] = await grok.dapi.tables.getTable(pulledRun.inputs[input.name]);
      return Promise.resolve();
    }));


    return pulledRun;
  }

  /**
   * Shallow saving already existing FuncCall with the same id. Should
   * be used only for updating run meta.
   * @param callToSave FuncCall to save
   */
  export async function updateRunMeta(call: DG.FuncCall) {
    await keepRunAuthor(call);
    const fc = await grok.dapi.functions.calls.allPackageVersions().save(call);
    grok.log.audit('@user edited metadata of @model model run', {
      'user': DG.User.current(),
      'model': fc.func,
    });
  }

  /**
   * Saved given FuncCall.
   * FuncCall is only stores references to actual dataframes. Thus, we should upload them separately
   * @param callToSave FuncCall to save
   * @param options audience: the group granted View on the uploaded dataframes (default: All users);
   * newId: assign the call a fresh id before saving, leaving an already saved run with the old id untouched;
   * metadata: the run's name, description and tags;
   * target: a space to link the run into (My runs otherwise), and the runs it replaces there
   * @returns Saved FuncCall
   */
  export async function saveRun(
    callToSave: DG.FuncCall,
    options?: {audience?: DG.Group, newId?: boolean, metadata?: RunMetadata, target?: SaveTarget, nested?: boolean},
  ) {
    let allGroup = groupsCache.get('All users');

    if (!allGroup) {
      allGroup = await grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
      groupsCache.set('All users', allGroup);
    }
    const audience = options?.audience ?? allGroup;

    if (options?.newId) callToSave.newId();
    if (options?.metadata) setRunMetadata(callToSave, options.metadata);
    callToSave.options[CREATED_BY_OPTION] = DG.User.current().id;
    callToSave.options[CREATED_AT_OPTION] = dayjs().toISOString();
    // a run saved on its own is not nested, even when it was loaded from a nested one
    if (options?.nested)
      callToSave.options[NESTED_ENTITY_OPTION] = 'true';
    else
      delete callToSave.options[NESTED_ENTITY_OPTION];
    const callCopy = deepCopy(callToSave);
    // a run opened from a space may have been deleted from its author's My runs; its copy is a new, live run
    delete callCopy.options['isDeleted'];
    if (options?.metadata) callCopy.friendlyName = callToSave.friendlyName;
    if (isIncomplete(callCopy)) callCopy.options['createdOn'] = dayjs().utc(true).unix();

    const fileInputs = wu(callCopy.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.FILE &&
        !!callCopy.inputs[input.name]);

    await Promise.all(fileInputs
      .map(async (input) => {
        // Compute2 stores file inputs as DG.FileInfo; compute1 stores a browser File.
        const value = callCopy.inputs[input.name];
        const bytes = value instanceof DG.FileInfo ?
          await value.readAsBytes() :
          new Uint8Array(await (value as File).arrayBuffer());
        const name = value instanceof DG.FileInfo ? (value.fileName ?? value.name) : (value as File).name;
        const filledFileInfo = DG.FileInfo.fromBytes(name, bytes);
        await grok.dapi.files.write(filledFileInfo);
        callCopy.inputs[input.name] = {id: filledFileInfo.id, name: filledFileInfo.name};

        return Promise.resolve();
      }),
    );

    const dfOutputs = wu(callCopy.outputParams.values() as DG.FuncCallParam[])
      .filter((output) =>
        output.property.propertyType === DG.TYPE.DATA_FRAME &&
        !!callCopy.outputs[output.name],
      );

    await Promise.all(dfOutputs
      .map(async (output) => {
        await grok.dapi.tables.uploadDataFrame(callCopy.outputs[output.name]);
        await grok.dapi.permissions.grant(callCopy.outputs[output.name].getTableInfo(), audience, false);
      }));

    const dfInputs = wu(callCopy.inputParams.values() as DG.FuncCallParam[])
      .filter((input) =>
        input.property.propertyType === DG.TYPE.DATA_FRAME &&
        !!callCopy.inputs[input.name],
      );

    await Promise.all(dfInputs
      .map(async (input) => {
        await grok.dapi.tables.uploadDataFrame(callCopy.inputs[input.name]);
        await grok.dapi.permissions.grant(callCopy.inputs[input.name].getTableInfo(), audience, false);
      }));

    const fc = await grok.dapi.functions.calls.allPackageVersions().save(callCopy);
    // await grok.dapi.permissions.grant(fc, audience, false);
    //   ^ requires core work: FuncCall is not a full entity — it has no per-run ACL,
    //     so only the uploaded dataframes are audience-scoped and the run itself
    //     stays loadable by id for everyone.
    // await grok.dapi.functions.calls.transferOwnership(fc, owner);
    //   ^ requires core work: ownership hand-off, so published/frozen copies stop
    //     being the author's to edit or soft-delete.
    grok.log.audit('@user saved run of the @model model', {
      'user': DG.User.current(),
      'model': fc.func,
    });
    if (options?.target?.spaceId) {
      try {
        await linkRun(fc.id, options.target.spaceId, options.target.replace);
      } catch (e: any) {
        // the run is saved, so the save still succeeds
        grok.shell.error(`Saved the run, but couldn't add it to the space: ${e?.message ?? e}`);
      }
    }
    // old views save without options and keep to their own history
    if (options && !options.nested)
      grok.events.fireCustomEvent(RUN_SAVED_EVENT, {spaceIds: options.target?.spaceId ? [options.target.spaceId] : []});
    return fc;
  }

  /** Links an existing run into a space, then removes the links of the runs it replaces there. */
  export async function linkRun(runId: string, spaceId: string, replace: DG.FuncCall[] = []) {
    const space = grok.dapi.spaces.id(spaceId);
    await space.addEntity(runId, true);
    for (const replaced of replace) {
      if (replaced.id !== runId)
        await space.removeEntity(replaced.id);
    }
  }

  /** Removes a run's link from a space. The run itself doesn't change. */
  export async function unlinkRun(runId: string, spaceId: string) {
    await grok.dapi.spaces.id(spaceId).removeEntity(runId);
  }

  /** All of the current user's runs of the models that the user hasn't deleted, wherever they are linked.
   * A model's runs include the runs saved from its custom view. */
  export async function listMyRuns(models: DG.Func | DG.Func[]): Promise<DG.FuncCall[]> {
    const funcs = Array.isArray(models) ? models : [models];
    if (funcs.length === 0)
      return [];
    const me = DG.User.current().id;
    const mine = `(options.${CREATED_BY_OPTION} = "${me}") or ` +
      `((session.user.id="${me}") and (options.${CREATED_BY_OPTION} = null))`;
    const names = [...new Set(funcs.map((func) => func.name))];
    const views = funcs.map((func) => `(options.${CUSTOM_VIEW_OPTION} = "${func.nqName}")`);
    const ofModels = `(func.name in (${names.map((name) => `"${name}"`).join(',')})) or ${views.join(' or ')}`;
    // every clause in its own parentheses, and the or-null clauses first: the server drops the first branch of a
    // nested or when an or-null clause comes before it unwrapped, and matches nothing when an options equality
    // comes before an or-null clause
    const runs = await grok.dapi.functions.calls.allPackageVersions()
      .filter(`(${mine}) and (${NOT_DELETED}) and (options.${NESTED_ENTITY_OPTION} = null) and (${ofModels})`)
      .include('session.user,options,func.package')
      .list();
    return runs.filter((run) => funcs.some((func) => isRunOf(func, run) &&
      (run.options?.[CUSTOM_VIEW_OPTION] === func.nqName || func.package?.name === run.func?.package?.name)));
  }

  /** Items owned by a space, and the runs linked into it. */
  export async function listSpace(spaceId: string): Promise<{items: DG.Entity[], runs: DG.FuncCall[]}> {
    const [items, runs] = await Promise.all([
      grok.dapi.spaces.id(spaceId).children.list(),
      listSpaceRuns(spaceId),
    ]);
    return {items, runs};
  }

  /** Runs linked into a space. */
  export async function listSpaceRuns(spaceId: string): Promise<DG.FuncCall[]> {
    // ofTypes changes the client it is called on, so each request gets its own
    return await grok.dapi.spaces.id(spaceId).children.ofTypes('FuncCall', true).list() as DG.FuncCall[];
  }

  /** Runs and the first other item in a space that already have this name, ignoring case. */
  export async function findNameConflict(name: string, spaceId: string): Promise<NameConflict> {
    return nameConflictIn(await listSpace(spaceId), name);
  }

  /** Whether a run is a run of the model, or was saved from the model's custom view; matched across packages and
   * versions by function name. */
  export const isRunOf = (func: DG.Func, run: DG.FuncCall) =>
    func.name === run.func?.name || func.nqName === run.options?.[CUSTOM_VIEW_OPTION];

  export const sameRunName = (a: string, b: string) => normalizeName(a) === normalizeName(b);

  /** A space item other than a run that users see and name; every space has a storage connection they don't. */
  export const isNamedSpaceItem = (item: DG.Entity) =>
    item.entityType !== 'FuncCall' && item.entityType !== 'DataConnection';

  /** The runs and the first other item of a space listing that have this name, ignoring case. */
  export function nameConflictIn(listing: {items: DG.Entity[], runs: DG.FuncCall[]}, name: string): NameConflict {
    const key = normalizeName(name);
    const other = listing.items.find((item) =>
      isNamedSpaceItem(item) && normalizeName(item.friendlyName || item.name || '') === key);
    return {runs: listing.runs.filter((run) => normalizeName(getRunName(run)) === key), other};
  }

  export async function deleteRun(callToDelete: DG.FuncCall) {
    await keepRunAuthor(callToDelete);
    callToDelete.options['isDeleted'] = true;
    await grok.dapi.functions.calls.allPackageVersions().save(callToDelete);
  }
}
