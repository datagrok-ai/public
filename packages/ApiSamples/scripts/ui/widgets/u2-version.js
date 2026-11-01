//api: DG.U2.version
// The version of the u2 engine the platform serves as DG.U2. A plugin built with u2 checks it when it loads: on an
// engine older than the one it needs, it fails at once with "u2 needs the u2 engine 1.2+ …: update the platform".
const [major, minor] = DG.U2.version.split('.').map(Number);
grok.shell.info(`u2 engine ${DG.U2.version}: ${major > 1 || major === 1 && minor >= 2 ? 'runs' : 'too old for'} u2`);
