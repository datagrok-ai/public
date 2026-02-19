/**
 * Progress indicator widget for the task bar.
 * @module widgets/progress
 */

import {toJs} from "../wrappers";
import {ProgressIndicator} from "../entities";
import {IDartApi} from "../api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class TaskBarProgressIndicator extends ProgressIndicator {
  static create(name?: string, options?: { cancelable?: boolean, pausable?: boolean, spinner?: boolean }): TaskBarProgressIndicator {
    return toJs(api.grok_TaskBarProgressIndicator_Create(name, options?.cancelable ?? false, options?.pausable ?? false, options?.spinner ?? false));
  }

  close(): any {
    return api.grok_TaskBarProgressIndicator_Close(this.dart);
  }
}
