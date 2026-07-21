/**
 * Node.js package loader — the headless equivalent of the browser's lazy package
 * bundle loading (xamgle ClientPackageFunc._loadPackage):
 *
 * 1. resolves the current published version via the REST API,
 * 2. fetches `dist/package.js` (a webpack `var`-library; its externals — DG, grok,
 *    ui, rxjs, dayjs, wu — are already globals after `startDatagrok`),
 * 3. evaluates it in this context and sets the `_package` fields,
 * 4. registers every package JS function with the local function registry
 *    (`grok.functions.register` → Node `RegisterFunc` interop), so
 *    `grok.functions.call('Pkg:fn')` runs it in-process with full param marshaling.
 *
 * UI-touching functions load fine and throw `DGNotSupportedError` only when called.
 */
import * as vm from 'node:vm';

const loaded: Map<string, any> = new Map();

export async function loadPackage(name: string): Promise<any> {
  const key = name.toLowerCase();
  if (loaded.has(key))
    return loaded.get(key);

  const g: any = globalThis;
  const grok = g.grok;
  if (!grok?.dapi?.root)
    throw new Error('loadPackage requires startDatagrok() to be called first');
  const root: string = grok.dapi.root;
  const headers = {Authorization: grok.dapi.token};

  // Resolve the published version the platform would serve: the caller's debug
  // version if one exists, else the current one (PackagesService.getCurrentPublishedPackages).
  const resp = await fetch(`${root}/packages?text=${encodeURIComponent(name)}`, {headers});
  if (!resp.ok)
    throw new Error(`Failed to list packages: HTTP ${resp.status}`);
  const packages: any[] = await resp.json();
  const pkg = packages.find((p) => p.name?.toLowerCase() === key);
  if (!pkg)
    throw new Error(`Package "${name}" is not published on ${root}`);
  const versions: any[] = pkg.publishedVersions ?? [];
  const login = (await grok.dapi.users.current()).login;
  const cur =
    versions.find((v) => v.debug && v.version === login) ??
    versions.find((v) => v.isCurrent) ??
    versions.find((v) => v.isLatest) ??
    versions[0];
  if (!cur)
    throw new Error(`Package "${name}" has no published versions`);

  const build = `${cur.version}/${cur.buildHash}/${cur.buildNumber}`;
  const webRoot = `${root}/packages/published/files/${pkg.name}/${build}/`;
  const jsUrl = `${webRoot}${cur.isWebpack === false ? '' : 'dist/'}package.js`;
  const codeResp = await fetch(jsUrl, {headers});
  if (!codeResp.ok)
    throw new Error(`Failed to fetch ${jsUrl}: HTTP ${codeResp.status}`);
  const code = await codeResp.text();

  // webpack var-library: top-level `var <name.toLowerCase()> = ...` lands on the global
  vm.runInThisContext(code, {filename: jsUrl});
  const module = g[key];
  if (module == null)
    throw new Error(`Bundle of "${name}" did not define the expected module global "${key}"`);

  if (module._package != null) {
    module._package.webRoot = webRoot;
    module._package.name = pkg.name;
    module._package.version = cur.version;
  }

  // Same registration the browser performs after a bundle loads, but with the real
  // output signature (the browser binds into ClientPackageFunc entities instead).
  const funcs: any[] = await grok.dapi.functions.filter(`package.shortName = "${pkg.name}"`).list();
  const registered: string[] = [];
  for (const f of funcs) {
    if ((f.options['file'] ?? 'package.js') !== 'package.js')
      continue;
    const impl = module[f.name];
    if (typeof impl !== 'function')
      continue;
    const ins = f.inputs.map((p: any) => `${p.propertyType} ${p.name}`).join(', ');
    const outs = f.outputs ?? [];
    const signature = outs.length > 1
      ? `({${outs.map((p: any) => `${p.propertyType} ${p.name}`).join(', ')}}) ${f.name}(${ins})`
      : `${outs[0]?.propertyType ?? 'bool'} ${f.name}(${ins})`;
    grok.functions.register({
      signature,
      isAsync: true,
      namespace: pkg.name,
      run: async (...args: any[]) => await impl(...args),
    });
    registered.push(f.name);
  }

  for (const f of funcs)
    if (f.options['tags']?.includes('init') || f.tags?.includes('init'))
      if (typeof module[f.name] === 'function')
        await module[f.name]();

  const result = {name: pkg.name, version: cur.version, webRoot, module, functions: registered};
  loaded.set(key, result);
  return result;
}
