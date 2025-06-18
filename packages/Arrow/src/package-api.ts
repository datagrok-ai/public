import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:Info', {});
  }

  export async function parquetInit(): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:ParquetInit', {});
  }

  export async function initPackage(): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:InitPackage', {});
  }

  //Converts DG.DataFrame to arrow
  export async function toFeather(table: DG.DataFrame, asStream: boolean): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:ToFeather', { table, asStream });
  }

  //Converts arrow ipc stream to DG.DataFrame
  export async function fromFeather(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:FromFeather', { bytes });
  }

  //Converts DG.DataFrame to parquet
  export async function toParquet(table: DG.DataFrame, compression: number): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:ToParquet', { table, compression });
  }

  //Converts binary data in parquet format to DG.DataFrame
  export async function fromParquet(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:FromParquet', { bytes });
  }

  export async function parquetFileHandler(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:ParquetFileHandler', { bytes });
  }

  export async function featherFileHandler(bytes: any): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:FeatherFileHandler', { bytes });
  }

  //Save as Parquet
  export async function saveAsParquet(): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:SaveAsParquet', {});
  }

  //Save as Feather
  export async function saveAsFeather(): Promise<any> {
    return await grok.functions.call('@datagrok/arrow:SaveAsFeather', {});
  }
}
