/// this file was generated automatically from grok_shared classes declarations
import { toDart } from "../wrappers";
let api = <any>window;

import {Entity} from '../entities'

export class DockerImage extends Entity {
  constructor(dart: any) { super(dart); };

  static fromJson(map: Map<any, any>): DockerImage { return new DockerImage(api.grok_DockerImage_fromJson(toDart(map)))};
  static STATUS_READY = 'ready';

  static STATUS_ERROR = 'error';

  static STATUS_BUILDING = 'building';

  static STATUS_PENDING_BUILD = 'pending build';

  static dbTableName = 'dockerfiles';

  get description(): string { return api.grok_DockerImage_Get_description(this.dart); };
  set description(x: string) {api.grok_DockerImage_Set_description(this.dart, toDart(x)); }
  get dockerfile(): string { return api.grok_DockerImage_Get_dockerfile(this.dart); };
  set dockerfile(x: string) {api.grok_DockerImage_Set_dockerfile(this.dart, toDart(x)); }
  get status(): string { return api.grok_DockerImage_Get_status(this.dart); };
  set status(x: string) {api.grok_DockerImage_Set_status(this.dart, toDart(x)); }
  get dockerName(): string { return api.grok_DockerImage_Get_dockerName(this.dart); };
  set dockerName(x: string) {api.grok_DockerImage_Set_dockerName(this.dart, toDart(x)); }
  get version(): string { return api.grok_DockerImage_Get_version(this.dart); };
  set version(x: string) {api.grok_DockerImage_Set_version(this.dart, toDart(x)); }
  get dockerfilePath(): string { return api.grok_DockerImage_Get_dockerfilePath(this.dart); };
  set dockerfilePath(x: string) {api.grok_DockerImage_Set_dockerfilePath(this.dart, toDart(x)); }
  get updatedBy(): string { return api.grok_DockerImage_Get_updatedBy(this.dart); };
  set updatedBy(x: string) {api.grok_DockerImage_Set_updatedBy(this.dart, toDart(x)); }
  get logs(): string { return api.grok_DockerImage_Get_logs(this.dart); };
  set logs(x: string) {api.grok_DockerImage_Set_logs(this.dart, toDart(x)); }
  get completed(): boolean { return api.grok_DockerImage_Get_completed(this.dart); };

  get iconStatus(): string { return api.grok_DockerImage_Get_iconStatus(this.dart); };

  get dockerFullName(): string { return api.grok_DockerImage_Get_dockerFullName(this.dart); };

}
