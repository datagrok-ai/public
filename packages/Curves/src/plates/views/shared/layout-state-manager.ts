import * as DG from 'datagrok-api/dg';
import {Plate} from '../../../plate/plate';
import {Subject, Observable} from 'rxjs';

/** Defines a layout pattern that can be applied to a plate. */
export interface LayoutPattern {
  name: string;
  // The apply function now takes a 'role' to assign to the wells.
  apply: (plate: Plate, role: string) => void;
}

/** Represents a single layout application action in the history. */
export interface LayoutAction {
  id: string;
  patternName: string;
  role: string; // The role assigned by this action
  timestamp: Date;
  // Store the state of the plate *before* this action was applied for undo purposes.
  previousState: DG.DataFrame;
}

export class LayoutStateManager {
  private _plate: Plate;
  private _history: LayoutAction[] = [];
  private _stateChange$ = new Subject<void>();

  constructor(initialPlate: Plate) {
    this._plate = initialPlate;
  }

  /** The main plate object being managed. */
  get plate(): Plate {
    return this._plate;
  }

  /** The history of all applied layouts. */
  get history(): readonly LayoutAction[] {
    return this._history;
  }

  /** Observable that fires whenever the plate state or history changes. */
  get onStateChange(): Observable<void> {
    return this._stateChange$.asObservable();
  }

  /** Applies a new layout pattern to the plate and records it in the history. */
  applyLayout(pattern: LayoutPattern, role: string): void {
    if (!role) {
      console.warn('Role cannot be empty.');
      return;
    }
    const action: LayoutAction = {
      id: `${Date.now()}-${Math.floor(Math.random() * 10000)}`,
      patternName: pattern.name,
      role: role,
      timestamp: new Date(),
      previousState: this._plate.data.clone(), // Save current state for undo
    };

    this._history.push(action);
    pattern.apply(this._plate, role); // Modify the plate using the specified role
    this._stateChange$.next();
  }

  /** Undoes a specific action from the history. */
  undoAction(actionId: string): void {
    const actionIndex = this._history.findIndex((a) => a.id === actionId);
    if (actionIndex === -1) return;

    // The state to restore is the one *before* the action was applied.
    const actionToUndo = this._history[actionIndex];
    this._plate.data = actionToUndo.previousState.clone();

    // Remove this action and all subsequent actions from the history.
    this._history.splice(actionIndex);

    this._stateChange$.next();
  }
}


