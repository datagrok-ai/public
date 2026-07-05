package grok_connect.table_mutation;

import java.sql.BatchUpdateException;
import java.sql.SQLException;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

/**
 * Pure unit tests for the SQLState -> structured error-code map (connector-writes WO-6). No Docker.
 */
class SqlStateMapperTest {
    @DisplayName("Known SQLStates map to structured codes")
    @Test
    public void knownStates() {
        Assertions.assertEquals("unique", SqlStateMapper.code("23505"));
        Assertions.assertEquals("fk", SqlStateMapper.code("23503"));
        Assertions.assertEquals("notnull", SqlStateMapper.code("23502"));
        Assertions.assertEquals("check", SqlStateMapper.code("23514"));
        Assertions.assertEquals("conflict", SqlStateMapper.code("40001"));
        Assertions.assertEquals("conflict", SqlStateMapper.code("40P01"));
    }

    @DisplayName("The 22xxx data-exception class maps to 'value'")
    @Test
    public void valueClass() {
        Assertions.assertEquals("value", SqlStateMapper.code("22001")); // string data right truncation
        Assertions.assertEquals("value", SqlStateMapper.code("22003")); // numeric value out of range
        Assertions.assertEquals("value", SqlStateMapper.code("22P02")); // invalid text representation
    }

    @DisplayName("Unknown / null / short SQLStates fall back to db-error")
    @Test
    public void fallback() {
        Assertions.assertEquals("db-error", SqlStateMapper.code(null));
        Assertions.assertEquals("db-error", SqlStateMapper.code(""));
        Assertions.assertEquals("db-error", SqlStateMapper.code("2"));
        Assertions.assertEquals("db-error", SqlStateMapper.code("42P01")); // undefined table
        Assertions.assertEquals("db-error", SqlStateMapper.code("08006")); // connection failure
        Assertions.assertEquals("db-error", SqlStateMapper.code("23000")); // integrity constraint (unmapped generic)
    }

    @DisplayName("A BatchUpdateException is unwrapped to its next-exception cause")
    @Test
    public void unwrapsBatchException() {
        SQLException cause = new SQLException("duplicate key", "23505");
        BatchUpdateException batch = new BatchUpdateException("batch failed", new int[] {1, -3});
        batch.setNextException(cause);
        Assertions.assertSame(cause, SqlStateMapper.unwrap(batch));
        Assertions.assertSame(cause, SqlStateMapper.unwrap(cause)); // plain exception is returned as-is
    }
}
