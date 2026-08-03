import * as rxjs from 'rxjs';
import {CardView} from './card_view';
import {Project} from '../entities/project';
import {Entity} from '../entities/entity';
import {__obs} from '../events';
import {toJs} from '../wrappers';
import {IDartApi} from '../api/grok_api.g';
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Hierarchical view of a Datagrok Space — the same widget the Browse tree
 *  opens when you click a Space (a {@link Project} where `isSpace == true`).
 *  Includes the search box, the "+" ribbon button to add a child space, and
 *  the card grid of children (subspaces, files, entities). */
export class SpaceView extends CardView {
  /** Builds the view bound to the given space. */
  static forProject(project: Project): SpaceView {
    return new SpaceView(api.grok_SpaceView_Create_ForProject(project.dart));
  }

  /** The space (Project) being shown. */
  get project(): Project { return toJs(api.grok_SpaceView_Get_Project(this.dart)); }

  /** Whether the in-view preview dock is shown when an item is selected.
   *  Set to `false` when you want to render the selected item's preview
   *  yourself by subscribing to {@link onCurrentObjectChanged}. */
  get showItemPreview(): boolean { return api.grok_SpaceView_Get_ShowItemPreview(this.dart); }
  set showItemPreview(s: boolean) { api.grok_SpaceView_Set_ShowItemPreview(this.dart, s); }

  /** Fires when the user selects a different item in the space's card grid. */
  get onCurrentObjectChanged(): rxjs.Observable<Entity> {
    return __obs<Entity>('grok-view-current-item-changed', this.dart);
  }
}
