export async function knownOpenBug(
  bugId: string,
  assertion: () => void | Promise<void>,
): Promise<void> {
  try {
    await assertion();
  } catch {
    console.log(`[KNOWN_BUG_REPRODUCED:${bugId}]`);
    return;
  }
  throw new Error(
    `[KNOWN_BUG_FIXED:${bugId}] assertion now passes — set related_bugs[].status: ` +
    `fixed + fixed_in and replace knownOpenBug with a plain hard expect`,
  );
}
