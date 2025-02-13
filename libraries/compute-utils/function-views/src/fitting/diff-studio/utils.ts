import {IVP} from '@datagrok/diff-studio-tools';

export function isWorkerApplicable(ivp: IVP): boolean {
  return (ivp.loop === null) && (ivp.updates === null) && (ivp.outputs === null);
}
