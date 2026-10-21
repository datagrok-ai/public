// The image starts the server as root and the server drops to uid 1001 here, before anything else runs: a process
// whose uid changed after exec is not dumpable, so the agents' tools (also uid 1001) cannot read its memory, environ
// or fds, which hold the provider key and every connected user's key.
if (process.getuid?.() === 0) {
  process.setgroups!([]);
  process.setgid!(2001);
  process.setuid!(1001);
}
