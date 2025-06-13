import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function exceptionScriptJulia(a: number): Promise<number> {
    return await grok.functions.call('@datagrok/dev-tools:ExceptionScriptJulia', { a });
  }

  export async function exceptionScriptOctave(a: number): Promise<number> {
    return await grok.functions.call('@datagrok/dev-tools:ExceptionScriptOctave', { a });
  }

  export async function exceptionScriptPython(a: number): Promise<number> {
    return await grok.functions.call('@datagrok/dev-tools:ExceptionScriptPython', { a });
  }
}

export namespace queries {
  export async function testHistory(packageName: string, category: string, test: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dev-tools:TestHistory', { packageName, category, test });
  }

  export async function categoryHistory(packageName: string, category: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dev-tools:CategoryHistory', { packageName, category });
  }

  export async function packageHistory(packageName: string): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/dev-tools:PackageHistory', { packageName });
  }
}

export namespace funcs {
  export async function renderDevPanel(ent: any): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:RenderDevPanel', { ent });
  }

  export async function makeInspectorPanel(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:MakeInspectorPanel', {});
  }

  //DevTools autostart function
  export async function autostartTools(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:AutostartTools', {});
  }

  //IconTool
  export async function iconTool(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:IconTool', {});
  }

  export async function testManager(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:TestManager', {});
  }

  export async function testDetectors(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:TestDetectors', {});
  }

  export async function testDetectorsStandard(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:TestDetectorsStandard', {});
  }

  export async function testFunctions(scope: any): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:TestFunctions', { scope });
  }

  export async function testFunction(): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:TestFunction', {});
  }

  export async function exceptionFunc(a: number): Promise<any> {
    return await grok.functions.call('@datagrok/dev-tools:ExceptionFunc', { a });
  }
}
