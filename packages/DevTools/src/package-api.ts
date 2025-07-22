import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function exceptionScriptJulia(a: number): Promise<number> {
    return await grok.functions.call('DevTools:ExceptionScriptJulia', { a });
  }

  export async function exceptionScriptOctave(a: number): Promise<number> {
    return await grok.functions.call('DevTools:ExceptionScriptOctave', { a });
  }

  export async function exceptionScriptPython(a: number): Promise<number> {
    return await grok.functions.call('DevTools:ExceptionScriptPython', { a });
  }
}

export namespace queries {
  export async function testHistory(packageName: string, category: string, test: string): Promise<DG.DataFrame> {
    return await grok.data.query('DevTools:TestHistory', { packageName, category, test });
  }

  export async function categoryHistory(packageName: string, category: string): Promise<DG.DataFrame> {
    return await grok.data.query('DevTools:CategoryHistory', { packageName, category });
  }

  export async function packageHistory(packageName: string): Promise<DG.DataFrame> {
    return await grok.data.query('DevTools:PackageHistory', { packageName });
  }
}

export namespace funcs {
  export async function renderDevPanel(ent: any): Promise<any> {
    return await grok.functions.call('DevTools:RenderDevPanel', { ent });
  }

  export async function makeInspectorPanel(): Promise<any> {
    return await grok.functions.call('DevTools:MakeInspectorPanel', {});
  }

  //DevTools autostart function
  export async function autostartTools(): Promise<any> {
    return await grok.functions.call('DevTools:AutostartTools', {});
  }

  //IconTool
  export async function iconTool(): Promise<any> {
    return await grok.functions.call('DevTools:IconTool', {});
  }

  export async function testManager(): Promise<any> {
    return await grok.functions.call('DevTools:TestManager', {});
  }

  export async function testDetectors(): Promise<any> {
    return await grok.functions.call('DevTools:TestDetectors', {});
  }

  export async function testDetectorsStandard(): Promise<any> {
    return await grok.functions.call('DevTools:TestDetectorsStandard', {});
  }

  export async function testFunctions(scope: any): Promise<any> {
    return await grok.functions.call('DevTools:TestFunctions', { scope });
  }

  export async function testFunction(): Promise<any> {
    return await grok.functions.call('DevTools:TestFunction', {});
  }

  export async function exceptionFunc(a: number): Promise<any> {
    return await grok.functions.call('DevTools:ExceptionFunc', { a });
  }
}
