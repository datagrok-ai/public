/// Docs: [Entity export / import](/docs/features/grok-tool/export-import/DESIGN.md)

/** Runs `work` over `items` with at most `concurrency` in flight, in the order given. */
export async function pool<T>(items: T[], concurrency: number, work: (item: T) => Promise<void>): Promise<void> {
  const queue = items.slice();
  const workers: Promise<void>[] = [];
  for (let i = 0; i < Math.min(concurrency, queue.length); i++)
    workers.push((async () => {
      while (queue.length)
        await work(queue.shift()!);
    })());
  await Promise.all(workers);
}
