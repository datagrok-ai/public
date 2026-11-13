/** The u2 engine's version, `DG.U2.version`. u2 refuses to build its controls on an engine older than the one it
 * needs, so a plugin built with a newer u2 than the platform carries says so instead of failing later on a missing
 * member, and still loads. */
export const version: string = '1.2.0';
