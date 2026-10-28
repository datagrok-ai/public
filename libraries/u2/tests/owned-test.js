/* `node:test`'s `test` for code that builds u2 controls: each test is their owner. The body runs under a root scope of
   the test's own, which goes when the test ends, so what it builds before its first `await` is owned, as in
   `u2test`; after an `await`, `testOwner().build(() => …)`. A test of what happens with no owner takes `node:test`'s. */
import {test as nodeTest} from 'node:test';
import {Scope} from 'datagrok-api/u2core';

let running;

export function test(name, options, body) {
  if (typeof options === 'function')
    [options, body] = [{}, options];
  if (typeof body !== 'function')
    return nodeTest(name, options, body);
  return nodeTest(name, options, async (t) => {
    const owner = running = Scope.root('test');
    try {
      return await owner.build(() => body(t));
    }
    finally {
      running = undefined;
      owner.dispose();
    }
  });
}

/** The running test's scope. */
export function testOwner() {
  if (running === undefined)
    throw new Error('testOwner() is called from a test body');
  return running;
}
