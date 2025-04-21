import { CardView } from './card_view';
import { toJs } from '../wrappers';
import {IDartApi} from '../api/grok_api.g';
import {TreeViewGroup} from '../widgets';
const api: IDartApi = <any>window;


/** File browser view. Contains the file tree, search panel and preview components. */
export class FilesView extends CardView {
  static create(options?: any): FilesView {
    return new FilesView(api.grok_FilesView_Create(options));
  }

  /** File tree visibility. */
  get showTree(): boolean { return toJs(api.grok_FilesView_Get_ShowTree(this.dart)); }
  set showTree(s: boolean) { api.grok_FilesView_Set_ShowTree(this.dart, s); }

  /** File preview visibility. */
  get showPreview(): boolean { return toJs(api.grok_FilesView_Get_ShowPreview(this.dart)); }
  set showPreview(s: boolean) { api.grok_FilesView_Set_ShowPreview(this.dart, s); }

  /** Search panel visibility. */
  get showSearch(): boolean { return toJs(api.grok_FilesView_Get_ShowSearch(this.dart)); }
  set showSearch(s: boolean) { api.grok_FilesView_Set_ShowSearch(this.dart, s); }

  /** Shows the tree component only. Similar to [FilesWidget]. */
  get showTreeOnly(): boolean { return this.showTree && !this.showPreview && !this.showSearch; }
  set showTreeOnly(s: boolean) {
    this.showTree = true;
    this.showPreview = !s;
    this.showSearch = !s;
  }

  /** Returns the files view tree. */
  get tree(): TreeViewGroup { return toJs(api.grok_FilesView_Get_Tree(this.dart)); }
}
