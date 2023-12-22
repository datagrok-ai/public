import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Base64} from 'js-base64';

import {DataProviderFunc} from '../utils/data-provider';

export type BiostructureData = {
  binary: boolean,
  ext: string,
  data: string | Uint8Array,
  options?: { name?: string, dataLabel?: string, },
}

export type BiostructureDataProviderFunc = DataProviderFunc<string, BiostructureData>;

/**
 * Viewer does not support a property of type extending String event completely compatible
 */
export class BiostructureDataJson extends String {
  protected constructor(value: string) {
    super(value);
  }

  static null = 'null';

  static fromData(src: BiostructureData): string {
    const dataStr = src.binary ? Base64.fromUint8Array(src.data as Uint8Array) : src.data as string;
    return JSON.stringify({data: dataStr, ext: src.ext, binary: src.binary, options: src.options});
  }

  static toData(src: string): BiostructureData {
    const data = JSON.parse(src as unknown as string);
    if (data) data.data = data.binary ? Base64.toUint8Array(data.data) : data.data;
    return data;
  }
}

