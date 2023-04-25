import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {IViewer} from './viewer';

export enum PositionHeight {
  Entropy = 'Entropy',
  full = '100%',
}

export class WebLogoPropsDefault {
  // -- Data --
  sequenceColumnName: string | null = null;
  startPositionName: string | null = null;
  endPositionName: string | null = null;
  skipEmptySequences: boolean = true;
  skipEmptyPositions: boolean = false;
  shrinkEmptyTail: boolean = true;
}

export type WebLogoProps = Required<WebLogoPropsDefault>;

export interface IWebLogoViewer extends WebLogoProps, IViewer {
  // Some methods to control the viewer
}

export const positionSeparator: string = ', ';
export const positionRe: RegExp = /(\d+)([A-Z]?)/;

export enum TAGS {
  positionNames = '.positionNames',
}
