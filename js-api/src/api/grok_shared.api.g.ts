/// this file was generated automatically from grok_shared classes declarations
let api = <any>window;

import {Entity} from '../entities'

export class DockerImage extends Entity {
  constructor(dart: any) { super(dart); };

  static fromJson(map: Map<any, any>): DockerImage { return new DockerImage(api.grok_DockerImage_fromJson(map))};
  static STATUS_READY = 'ready';

  static STATUS_USED = 'used';

  static STATUS_ERROR = 'error';

  static STATUS_BUILDING = 'building';

  static STATUS_PENDING_BUILD = 'pending build';

  static dbTableName = 'dockerfiles';

  get description(): string { return api.grok_DockerImage_Get_description(this.dart); };
  set description(x: string) {api.grok_DockerImage_Set_description(this.dart, x); }
  get dockerfile(): string { return api.grok_DockerImage_Get_dockerfile(this.dart); };
  set dockerfile(x: string) {api.grok_DockerImage_Set_dockerfile(this.dart, x); }
  get status(): string { return api.grok_DockerImage_Get_status(this.dart); };
  set status(x: string) {api.grok_DockerImage_Set_status(this.dart, x); }
  get dockerName(): string { return api.grok_DockerImage_Get_dockerName(this.dart); };
  set dockerName(x: string) {api.grok_DockerImage_Set_dockerName(this.dart, x); }
  get version(): string { return api.grok_DockerImage_Get_version(this.dart); };
  set version(x: string) {api.grok_DockerImage_Set_version(this.dart, x); }
  get dockerfilePath(): string { return api.grok_DockerImage_Get_dockerfilePath(this.dart); };
  set dockerfilePath(x: string) {api.grok_DockerImage_Set_dockerfilePath(this.dart, x); }
  get className(): string { return api.grok_DockerImage_Get_className(this.dart); };

  get completed(): boolean { return api.grok_DockerImage_Get_completed(this.dart); };

  get iconStatus(): string { return api.grok_DockerImage_Get_iconStatus(this.dart); };

  get dockerFullName(): string { return api.grok_DockerImage_Get_dockerFullName(this.dart); };

}
