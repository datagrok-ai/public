import * as DG from 'datagrok-api/dg';


export function isApp(entity: DG.Entity): entity is DG.Func {
  return entity instanceof DG.Func && entity.options['role'] === DG.FUNC_TYPES.APP;
}

/** A model — a Func with the `model` role/tag. Launched via `prepare().edit()`, not run directly. */
export function isModel(entity: DG.Entity): entity is DG.Func {
  return entity instanceof DG.Func &&
    ((((entity.options['role'] as string) ?? '').split(',').includes('model')) || entity.hasTag('model'));
}

export function isRunnable(entity: DG.Entity): entity is DG.Func {
  return entity instanceof DG.Func && !isApp(entity) && !isModel(entity);
}

/** Returns true if the entity is relevant for Spotlight (Recent / Shared with me / Workspace). */
export function isSpotlightEntity(ent: DG.Entity): boolean {
  if (!ent || !ent.friendlyName)
    return false;
  if (ent instanceof DG.FuncCall || ent instanceof DG.Group || ent instanceof DG.User || ent instanceof DG.Package ||
    ent instanceof DG.UserReport || ent.entityType === 'UserReport' || ent instanceof DG.TableInfo ||
    (ent instanceof DG.Func && !(ent instanceof DG.Script || ent instanceof DG.DataQuery || ent instanceof DG.DataJob || isApp(ent))) ||
    ent instanceof DG.ViewInfo || ent instanceof DG.ViewLayout || ent instanceof DG.DataConnection ||
    (ent instanceof DG.Project && (ent.isPackage || (!ent.isDashboard && !ent.isSpace))) ||
    //@ts-ignore
      (ent.hasOwnProperty('npmScope') && ent['npmScope'] == 'datagrok'))
    return false;
  return true;
}
