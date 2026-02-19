/**
 * User, UserSession, and Group classes.
 * These are kept together due to circular dependencies (User.group → Group, Group.user → User).
 * @module entities/user
 */

import {USER_STATUS} from "../const";
import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import dayjs from "dayjs";

declare var grok: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

// Forward declarations for circular references
type DataConnection = any;
type Project = any;

/**
 * Represents a user of the Datagrok platform.
 * @extends Entity
 * */
export class User extends Entity {
  /** @constructs User*/
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new user.
   * Note that it's just a client object, it won't be saved in the database. */
  static create(): User { return new User(api.grok_User_From_Id(null)); };

  static fromId(id: string): User { return new User(api.grok_User_From_Id(id)); }

  /** Returns current user. */
  static current(): User { return new User(api.grok_User()); }

  /** First name */
  get firstName(): string { return api.grok_User_Get_FirstName(this.dart); }
  set firstName(name: string) {api.grok_User_Set_FirstName(this.dart, name);}

  /** Last name */
  get lastName(): string { return api.grok_User_Get_LastName(this.dart); }
  set lastName(name: string) {api.grok_User_Set_LastName(this.dart, name);}

  get status(): string & USER_STATUS { return api.grok_User_Get_Status(this.dart); }
  set status(name: string & USER_STATUS) {api.grok_User_Set_Status(this.dart, name);}

  /** Email */
  get email(): string | null { return api.grok_User_Get_Email(this.dart); }
  set email(email: string | null) {api.grok_User_Set_Email(this.dart, email);}

  /** Picture URL */
  get picture(): string | object { return api.grok_User_Get_Picture(this.dart); }

  /** User home project */
  get project(): Project { return toJs(api.grok_User_Get_Project(this.dart)); }

  /** User home folder connection */
  get home(): DataConnection { return toJs(api.grok_User_Get_Storage(this.dart)); }

  /** Login */
  get login(): string { return api.grok_User_Get_Login(this.dart); }
  set login(login: string) { api.grok_User_Set_Login(this.dart, login); }

  toMarkup(): string { return api.grok_User_ToMarkup(this.dart); }

  /** Security Group */
  get group(): Group { return toJs(api.grok_User_Get_Group(this.dart)); }

  /** Date when user joined */
  get joined(): dayjs.Dayjs { return dayjs(api.grok_User_Get_Joined(this.dart)); }

  static get defaultUsersIds() {
    return {
      "Test": "ca1e672e-e3be-40e0-b79b-d2c68e68d380",
      "Admin": "878c42b0-9a50-11e6-c537-6bf8e9ab02ee",
      "System": "3e32c5fa-ac9c-4d39-8b4b-4db3e576b3c3",
    } as const;
  }

  static get test(): User { return new User(api.grok_User_Test()); }

  static get admin(): User { return new User(api.grok_User_Admin()); }

  static get system(): User { return new User(api.grok_User_System()); }
}


/**
 * Represents a user session in the Datagrok platform.
 * @extends Entity
 * */
export class UserSession extends Entity {
  /** @constructs UserSession*/
  constructor(dart: any) {
    super(dart);
  }

  /** Entity ID (GUID) */
  get id(): string { return api.grok_Entity_Get_Id(this.dart); }
  set id(x: string) { api.grok_Entity_Set_Id(this.dart, x); }

  /** */
  get type(): 'internal' | 'guest' { return api.grok_UserSession_Get_Type(this.dart); }

  /** External Token */
  get externalToken(): string { return api.grok_UserSession_Get_ExternalToken(this.dart); }

  /** User */
  get user(): User { return toJs(api.grok_UserSession_Get_User(this.dart)); }
}


/** @extends Entity
 * Represents a User Group
 * */
export class Group extends Entity {
  /** @constructs Group */
  constructor(dart: any) {
    super(dart);
  }

  static create(name: string): Group { return new Group(api.grok_Group(name)); }

  /** Adds a member to the group
   * @param {Group} m */
  addMember(m: Group): void { api.grok_Group_Add_Member(this.dart, m.dart, false); }

  /** Adds an admin member to the group
   * @param {Group} m */
  addAdminMember(m: Group): void { api.grok_Group_Add_Member(this.dart, m.dart, true); }

  /** Removes a member from the group
   * @param {Group} m */
  removeMember(m: Group): void { api.grok_Group_Remove_Member(this.dart, m.dart); }

  /** Adds the group to another one
   * @param {Group} m */
  includeTo(m: Group): void { api.grok_Group_Add_Membership(this.dart, m.dart, false); }

  /** Adds the group to another one as an admin
   * @param {Group} m */
  includeAdminTo(m: Group): void { api.grok_Group_Add_Membership(this.dart, m.dart, true); }

  /** Removes membership from another group
   * @param {Group} m */
  excludeFrom(m: Group): void { api.grok_Group_Remove_Membership(this.dart, m.dart); }

  /** Returns list of groups that belong to group, with no admin permissions
   * @type {Array<Group>} */
  get members(): Group[] { return toJs(api.grok_Group_Get_Members(this.dart, false)); }

  /** Returns list of groups that belong to group, with admin permissions
   * @type {Array<Group>} */
  get adminMembers(): Group[] { return toJs(api.grok_Group_Get_Members(this.dart, true)); }

  /** Returns list of groups that group belongs to, with no admin permissions
   * @type {Array<Group>} */
  get memberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.dart, false)); }

  /** Returns list of groups that group belongs to, with admin permissions
   * @type {list<Group>} */
  get adminMemberships(): Group[] { return toJs(api.grok_Group_Get_Memberships(this.dart, true)); }

  /** Personal user group */
  get personal(): boolean { return api.grok_Group_Get_Personal(this.dart); }
  set personal(e: boolean) { api.grok_Group_Set_Personal(this.dart, e); }

  /** Hidden group */
  get hidden(): boolean { return api.grok_Group_Get_Hidden(this.dart); }
  set hidden(e: boolean) { api.grok_Group_Set_Hidden(this.dart, e); }

  /** Returns associated user.
   * Returns `null` if the group is not {@link personal}.
   * See [Groups](https://datagrok.ai/help/govern/access-control/users-and-groups#groups)
   */
  get user(): User { return new User(api.grok_Group_Get_User(this.dart)); }

  static get defaultGroupsIds() {
    return {
      "All users": "a4b45840-9a50-11e6-9cc9-8546b8bf62e6",
      "Developers": "ba9cd191-9a50-11e6-9cc9-910bf827f0ab",
      "Need to create": "00000000-0000-0000-0000-000000000000",
      "Test": "ca1e672e-e3be-40e0-b79b-8546b8bf62e6",
      "Admin": "a4b45840-9a50-11e6-c537-6bf8e9ab02ee",
      "System": "a4b45840-ac9c-4d39-8b4b-4db3e576b3c3",
      "Administrators": "1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5",
    } as const;
  }

  static get allUsers(): Group { return new Group(api.grok_Group_AllUsers()); }

  static get developers(): Group { return new Group(api.grok_Group_Developers()); }

  static get needToCreate(): Group { return new Group(api.grok_Group_NeedToCreate()); }

  static get test(): Group { return new Group(api.grok_Group_Test()); }

  static get admin(): Group { return new Group(api.grok_Group_Admin()); }

  static get system(): Group { return new Group(api.grok_Group_System()); }

  static get administrators(): Group { return new Group(api.grok_Group_Administrators()); }
}
