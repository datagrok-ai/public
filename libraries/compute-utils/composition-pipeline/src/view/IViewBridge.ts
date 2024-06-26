import {BehaviorSubject, Subject} from 'rxjs';
import {HooksRunner} from './HooksRunner';
import {ActionTriggered, Actions, ViewConfig, ViewConfigChanges, ViewValidationResult} from './ViewCommunication';

export interface IViewBridge {
  actionsInput: BehaviorSubject<Actions | undefined>;
  actionsEvents: Subject<ActionTriggered>;

  configInput: BehaviorSubject<ViewConfig | undefined>;
  configChanges: Subject<ViewConfigChanges>;

  validationInput: BehaviorSubject<ViewValidationResult | undefined>;

  blockRunInput: BehaviorSubject<boolean>;
}
