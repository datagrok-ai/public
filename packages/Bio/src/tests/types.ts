import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type DfReaderFunc = () => Promise<DG.DataFrame>;

export type ConverterFunc = (srcCol: DG.Column) => DG.Column;
