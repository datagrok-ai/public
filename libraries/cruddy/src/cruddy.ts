import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


export abstract class Entity {
  abstract get type(): string;
  abstract get properties(): DG.Property[];
}


export abstract class EntityCrud<TEntity extends Entity> {

  /** Creates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract create(entity: TEntity): Promise<TEntity>;

  /** Reads the entities, according to the specified filter. */
  abstract read(filter: {[name: string]: string}): Promise<DG.DataFrame>;

  /** Updates the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract update(entity: TEntity): Promise<TEntity>;

  /** Deletes the entity in the database.
   * Returns the entity on success, or throws an exception otherwise. */
  abstract delete(entity: TEntity): Promise<TEntity>;
}


export class DbQueryEntityCrud<TEntity extends Entity> extends EntityCrud<TEntity> {
  create(entity: TEntity): Promise<TEntity> {
    //return grok.data.query()
  }
}


export abstract class EntityView<TEntity extends Entity> {

}


export class App {

}