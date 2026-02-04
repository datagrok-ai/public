/**
 * UserReport, UserReportsRule, and UserNotification classes.
 * @module entities/reports
 */

import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import {User} from "./user";
import dayjs from "dayjs";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class UserReport extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  get isResolved(): boolean {
    return api.grok_UserReport_IsResolved(this.dart);
  }

  get jiraTicket(): string {
    return api.grok_UserReport_JiraTicket(this.dart);
  }

  get assignee(): User {
    return toJs(api.grok_UserReport_Assignee(this.dart));
  }

  get reporter(): User {
    return toJs(api.grok_UserReport_Reporter(this.dart));
  }

  get description(): string {
    return toJs(api.grok_UserReport_Description(this.dart));
  }

  get createdOn(): dayjs.Dayjs {
    return dayjs(api.grok_UserReport_CreatedOn(this.dart));
  }
}


export class UserReportsRule extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  static async showAddDialog(): Promise<void> {
    await api.grok_ReportsRule_Add_Dialog();
  }
}

export class UserNotification {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  get user(): User {
    return toJs(api.grok_UserNotification_User(this.dart));
  }

  get name(): string {
    return toJs(api.grok_UserNotification_Name(this.dart));
  }

  get friendlyName(): string {
    return toJs(api.grok_UserNotification_FriendlyName(this.dart));
  }

  get text(): string {
    return toJs(api.grok_UserNotification_Text(this.dart));
  }

  get data(): string {
    return toJs(api.grok_UserNotification_Data(this.dart));
  }

  get sender(): string {
    return toJs(api.grok_UserNotification_Sender(this.dart));
  }

  get createdAt(): dayjs.Dayjs {
    return dayjs(api.grok_UserNotification_CreatedAt(this.dart));
  }

  get isRead(): boolean {
    return api.grok_UserNotification_IsRead(this.dart);
  }

  get readAt(): dayjs.Dayjs {
    return dayjs(api.grok_UserNotification_ReadAt(this.dart));
  }
}
