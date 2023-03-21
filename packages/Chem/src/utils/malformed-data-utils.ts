import BitArray from "@datagrok-libraries/utils/src/bit-array";
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export function malformedDataWarning(fingerprintCol: BitArray[], dataframe: DG.DataFrame): (number | null)[] {
    const malformedData = fingerprintCol.map((it: BitArray, idx: number) => it.allFalse ? idx : null).filter((it) => it !== null);
      if (malformedData.length) {
        const message = ui.divText(`${malformedData.length} molecules with indexes: ${malformedData.map((it: number | null) => it! + 1)
          .join(',')} are empty or possibly malformed and are not included in analysis`);
        const selectRowsButton = ui.button('Select', ()=> {
          for (const i of malformedData)
          dataframe.selection.set(i!, true);
        })
        grok.shell.warning(ui.div([message, selectRowsButton]));
      }
    return malformedData;
}
