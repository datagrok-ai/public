package grok_connect.table_mutation;

import serialization.DataFrame;

/**
 * Consumes a streamed typed-DataFrame row payload for one bulk mutation (java-d42-reader WO-4).
 * Implementations own their own flushing; {@link MutationManager} owns the transaction (commit after
 * {@link #finish()}, rollback in {@link #abort()}). All calls happen on the single WebSocket message
 * thread — the wire protocol is strict ping-pong, so implementations need no internal synchronization.
 */
public interface BulkLoader {
    /** Consumes one decoded d42 chunk (already read off the copied Jetty buffer). */
    void feed(DataFrame chunk) throws Exception;

    /**
     * Consumes one CSV byte chunk (legacy transport, {@code payloadFormat == "csv"}). Throwing by
     * default — only {@link BatchInsertBulkLoader} still supports it; deleted in WO-7.
     */
    default void feedCsv(byte[] csvChunk) throws Exception {
        throw new UnsupportedOperationException("CSV payload is not supported by this loader; send payloadFormat='d42'");
    }

    /** Flushes the remainder and returns the affected-row count. */
    MutationResult finish() throws Exception;

    /** Releases the loader without committing (the caller rolls back). Must not throw. */
    void abort();
}
