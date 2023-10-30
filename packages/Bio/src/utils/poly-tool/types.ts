import {TRANSFORMATION_TYPE, CYCLIZATION_TYPE} from './const';

export type MetaData = {
  leftTerminal: string,
  rightTerminal: string,
  transformationType: TRANSFORMATION_TYPE,
  cyclizationType: CYCLIZATION_TYPE,
}

export type ConnectionData = {
  monomerPosition: number,
  attachmentPoint: number,
}

