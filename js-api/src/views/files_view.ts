import { View } from './view';
import { toJs } from '../wrappers';
let api = <any>window;


/** File browser view. Contains the file tree, search panel and preview components. */
export class FilesView extends View {
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

  refresh() { api.grok_FilesView_Refresh(this.dart); }
}
