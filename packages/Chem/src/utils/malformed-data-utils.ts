import BitArray from "@datagrok-libraries/utils/src/bit-array";
import * as grok from 'datagrok-api/grok';

export function malformedDataWarning(fingerprintCol: BitArray[]): (number | null)[] {
    const malformedData = fingerprintCol.map((it: BitArray, idx: number) => it.allFalse ? idx : null).filter((it) => it !== null);
      if (malformedData.length) {
        grok.shell.warning(`${malformedData.length} molecules with indexes: ${malformedData.map((it: number | null) => it! + 1)
          .join(',')} are empty or possibly malformed and are not included in analysis`);
      }
    return malformedData;
}