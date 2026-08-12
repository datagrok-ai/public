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

    /** Flushes the remainder and returns the affected-row count. */
    MutationResult finish() throws Exception;

    /** Releases the loader without committing (the caller rolls back). Must not throw. */
    void abort();
}
