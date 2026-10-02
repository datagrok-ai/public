import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';

export class ActivityCliffsTutorial extends Tutorial {
  get name() {
    return 'Activity Cliffs';
  }

  get description() {
    return 'Detects pairs of molecules ' +
      'with similar structures but different activity.';
  }

  get steps() {return 12;}
}
