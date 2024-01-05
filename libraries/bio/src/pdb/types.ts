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

/** Return {@link BiostructureData} (jsonified to {@link string}) by identifier value of type {@link string}. */
export type BiostructureDataProviderFunc = DataProviderFunc<string, string>;

export namespace BiostructureDataJson {
  /* avoid static methods */

  export const empty: string = 'null';

  /** {@link DG.Viewer} does not support a property of type extending String even completely compatible */
  export function fromData(src: BiostructureData): string {
    const dataStr = src.binary ? Base64.fromUint8Array(src.data as Uint8Array) : src.data as string;
    return JSON.stringify({data: dataStr, ext: src.ext, binary: src.binary, options: src.options});
  }

  export function toData(src: string): BiostructureData {
    const data = JSON.parse(src as unknown as string);
    if (data) data.data = data.binary ? Base64.toUint8Array(data.data) : data.data;
    return data;
  }
}
