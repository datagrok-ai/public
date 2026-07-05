package grok_connect.table_mutation;

/**
 * Consumes a streamed CSV row payload for one bulk mutation (connector-writes WO-5). Implementations
 * own their own flushing; {@link MutationManager} owns the transaction (commit after {@link #finish()},
 * rollback in {@link #abort()}). All calls happen on the single WebSocket message thread — the wire
 * protocol is strict ping-pong, so implementations need no internal synchronization.
 */
public interface BulkLoader {
    /** Consumes one CSV byte chunk. The array is owned by the loader (already copied off the Jetty buffer). */
    void feed(byte[] csvChunk) throws Exception;

    /** Flushes the remainder and returns the affected-row count. */
    MutationResult finish() throws Exception;

    /** Releases the loader without committing (the caller rolls back). Must not throw. */
    void abort();
}
