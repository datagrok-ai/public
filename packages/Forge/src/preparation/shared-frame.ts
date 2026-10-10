import * as DG from 'datagrok-api/dg';

/** A frame of [columns] themselves, without a copy; give it back with {@link releaseFrame} after use. */
export function sharedFrame(columns: DG.Column[]): DG.DataFrame {
  return DG.DataFrame.fromColumns(columns);
}

/** Removes every column from [frame]: a shared column stops having the frame as a parent (and sending it events). */
export function releaseFrame(frame: DG.DataFrame): void {
  for (const name of frame.columns.names())
    frame.columns.remove(name, false);
}

/** Runs [call] on a {@link sharedFrame} of [columns], given back when the call ends, whatever happens. */
export async function onFrame<T>(columns: DG.Column[], call: (frame: DG.DataFrame) => Promise<T>): Promise<T> {
  const frame = sharedFrame(columns);
  try {
    return await call(frame);
  } finally {
    releaseFrame(frame);
  }
}
