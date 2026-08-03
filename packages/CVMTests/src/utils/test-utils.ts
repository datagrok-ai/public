// Tricky strings reused for escaping-fidelity checks across script and celery tests.
export const escapingTestStrings = ['\t\n\t\tsdfdsf\t', ' sdfds \\\'\"""', ' \n ', '\'\""\'', '\n and \\n',
  String.raw`CO\C1=C(C(=C(C=C1)/C=N\N=C(N)N)Cl)OC`, '"', '\'', '\n', '\t', '\\', '\\n', '\\r', '\\t'];

export function randomString(length: number, chars: string) {
  let result = '';
  for (let i = length; i > 0; --i) result += chars[Math.round(Math.random() * (chars.length - 1))];
  return result;
}

export function isEqualBytes(bytes1: Uint8Array, bytes2: Uint8Array): boolean {
  if (bytes1.length !== bytes2.length)
    return false;

  for (let i = 0; i < bytes1.length; i++) {
    if (bytes1[i] !== bytes2[i])
      return false;
  }

  return true;
}
