/** Decompresses a gzip stream (a column section after the 0xCAFEBEEF marker and two lengths). */
export async function gunzip(data: Uint8Array, expectedLength?: number): Promise<Uint8Array> {
  const stream = new Blob([data.slice()]).stream().pipeThrough(new DecompressionStream('gzip'));
  const out = new Uint8Array(await new Response(stream).arrayBuffer());
  if (expectedLength !== undefined && out.length !== expectedLength)
    throw new Error(`Decompressed ${out.length} bytes, expected ${expectedLength}`);
  return out;
}
