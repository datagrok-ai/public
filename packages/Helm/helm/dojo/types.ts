import type {DojoType, DojoxType} from '@datagrok-libraries/js-draw-lite/src/types/dojo';

export type DojoWindowType = {
  dojo$: {
    /** dojo source file postfix: '' for compressed, '.uncompressed.js' for uncompressed */
    uncompressed: string;
    ctx: { <T>(id: string): T };
    initPromise: Promise<void>;
  },
  dojo: DojoType;
  dojox: DojoxType;
}
