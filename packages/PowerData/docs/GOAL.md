> **Current direction (see `STATUS.md`):** Shipping the pure-TypeScript
> streaming reader for now. WASM / Arrow C++ path is parked until we have a
> measured bottleneck that justifies the build cost. The text below is the
> original long-term vision.

Develop PowerData, a Datagrok plugin that would provide ultra-fast support for importing large (1-2GB) CSV files. Use Apache Arrow's CSV reader (C++, WASM-compiled) - just the CSV part. Key points: streaming everywhere (we never want to read the whole file in memory - instead, we should stream-chunk it to the wasm webworker, where we would carefully construct the dataframe), dictionary encoding.

With progress updates, and cancellation mechanism.
Initially, it would be enough to expose a function to the main menu that would let user pick a CSV file from disk, open it and add to the shell (grok.sheel.addTableView). Later on we will convert it to the proper importer.

Use existing instructions for working with mixed wasm/TypeScript plugins.

---

Key design decisions:
Streaming vs. bulk. Arrow's CSV reader supports streaming (read N rows at a time). For 2GB files you want this — it lets you show progressive loading in the UI and keeps peak WASM memory well below 2x file size. You'd emit d42 chunks as they're ready.
WASM memory pool. You're right that the worker's WASM instance gets its own linear memory, so a 2GB parse won't blow up the main thread's heap. But watch out — WASM linear memory can only grow, never shrink. You'll want to terminate and recreate the worker after large imports to reclaim that memory.
Dictionary encoding. Arrow already does dictionary encoding for string columns if you configure it. If your d42 format uses a compatible dictionary scheme, you can literally copy the dictionary and index arrays out of Arrow's buffers. This is where the real performance win is — you avoid re-scanning strings on the main thread.

---

Test files:
C:\data\4gb-smiles.csv
C:\data\2gb-smiles.csv