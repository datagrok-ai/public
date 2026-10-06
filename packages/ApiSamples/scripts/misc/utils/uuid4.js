// A random v4 UUID. Unlike crypto.randomUUID() this also works outside a secure context,
// so it is safe on a stand served over plain HTTP.

grok.shell.info(DG.Utils.uuid4()); // 3f1c2a9e-5b74-4d0a-9c31-7e6f8a2b1d45
