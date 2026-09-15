/* A v4 uuid without `crypto.randomUUID`, which is secure-context only: a stand served over plain
   HTTP on a hostname has no `randomUUID`, and the first draft row would throw. `getRandomValues`
   is available everywhere, and these ids — draft keys, transaction stamps — carry no secret. */

const HEX: string[] = [];
for (let i = 0; i < 256; i++)
  HEX.push(i.toString(16).padStart(2, '0'));

export function uuid4(): string {
  const b = crypto.getRandomValues(new Uint8Array(16));
  b[6] = (b[6] & 0x0f) | 0x40;
  b[8] = (b[8] & 0x3f) | 0x80;
  const h = (from: number, to: number): string => {
    let s = '';
    for (let i = from; i < to; i++)
      s += HEX[b[i]];
    return s;
  };
  return `${h(0, 4)}-${h(4, 6)}-${h(6, 8)}-${h(8, 10)}-${h(10, 16)}`;
}
