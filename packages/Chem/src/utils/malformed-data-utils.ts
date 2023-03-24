import BitArray from "@datagrok-libraries/utils/src/bit-array";
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export function malformedDataWarning(fingerprintCol: (BitArray | null)[], dataframe: DG.DataFrame): number[] {
  const malformedData: number[] = [];
  for (let i = 0; i < fingerprintCol.length; i++) {
    if (!fingerprintCol[i])
      malformedData.push(i);
  }
  if (malformedData.length) {
    const message = ui.divText(`${malformedData.length} molecules are possibly malformed and are not included in analysis`);
    const selectRowsButton = ui.button('Select', () => {
      for (const i of malformedData)
        dataframe.selection.set(i!, true);
    })
    grok.shell.warning(ui.div([message, selectRowsButton]));
  }
  return malformedData as number[];
}

export function setEmptyBitArraysForMalformed(fingerprintCol: (BitArray | null)[]): void {
  for (let i = 0; i < fingerprintCol.length; i++) {
    if (!fingerprintCol[i])
      fingerprintCol[i] = BitArray.fromBytes(new Uint8Array());
  }
}
