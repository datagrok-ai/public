/**
 * Logging related classes: LogEventType, LogEvent, LogEventParameter, LogEventParameterValue.
 * These are kept together due to circular references.
 * @module entities/logging
 */

import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import {UserSession} from "./user";
import dayjs from "dayjs";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class LogEventType extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Friendly name of the event type */
  get name(): string { return api.grok_LogEventType_Get_Name(this.dart); }

  get comment(): string { return api.grok_LogEventType_Get_Comment(this.dart); }
  set comment(comment) { api.grok_LogEventType_Set_Comment(this.dart, comment); }

  get isError(): boolean { return api.grok_LogEventType_Get_IsError(this.dart); }
  set isError(isError) { api.grok_LogEventType_Set_IsError(this.dart, isError); }

}

export class LogEvent extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Description of the event */
  get description(): string { return api.grok_LogEvent_Get_Description(this.dart); }

  /** Friendly name of the event */
  get name(): string { return api.grok_LogEvent_Get_Name(this.dart); }

  /** Session id of the event */
  get session(): UserSession | string { return toJs(api.grok_LogEvent_Get_Session(this.dart)); }

  /** Parameters of the event
   * @type {Array<LogEventParameterValue>} */
  get parameters(): LogEventParameterValue[] { return toJs(api.grok_LogEvent_Get_Parameters(this.dart)); }

  /** Type of the event
   * @type {LogEventType} */
  get eventType(): LogEventType { return toJs(api.grok_LogEvent_Get_Type(this.dart)); }

  get eventTime(): dayjs.Dayjs { return dayjs(api.grok_LogEvent_Get_EventTime(this.dart)); }
}

export class LogEventParameter extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Name of the parameter */
  get name(): string { return api.grok_LogEventParameter_Get_Name(this.dart); }

  /** Type of the parameter */
  get type(): string { return api.grok_LogEventParameter_Get_Type(this.dart); }
}

export class LogEventParameterValue extends Entity {
  constructor(dart: any) {
    super(dart);
  }

  /** Event of the parameter value
   * @type {LogEvent} */
  get event(): LogEvent {
    return toJs(api.grok_LogEventParameterValue_Get_Event(this.dart));
  }

  /** Parameter of the parameter value
   * @type {LogEventParameter} */
  get parameter(): LogEventParameter {
    return toJs(api.grok_LogEventParameterValue_Get_Parameter(this.dart));
  }

  /** Parameter value
   * @type {string} */
  get value(): string {
    return api.grok_LogEventParameterValue_Get_Value(this.dart);
  }
}
