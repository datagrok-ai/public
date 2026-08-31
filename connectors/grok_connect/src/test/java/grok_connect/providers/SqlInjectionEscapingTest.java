package grok_connect.providers;

import grok_connect.connectors_info.FuncParam;
import grok_connect.table_mutation.MutationValidationException;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import serialization.Types;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;

public class SqlInjectionEscapingTest {
    private final PostgresDataProvider postgres = new PostgresDataProvider();
    private final MsSqlDataProvider mssql = new MsSqlDataProvider();

    @Test
    @DisplayName("escapeSqlString doubles single quotes")
    public void escapeSqlStringDoublesQuotes() {
        assertEquals("O''Brien", JdbcDataProvider.escapeSqlString("O'Brien"));
        assertEquals("'' OR ''1''=''1", JdbcDataProvider.escapeSqlString("' OR '1'='1"));
        assertEquals("plain", JdbcDataProvider.escapeSqlString("plain"));
    }

    @Test
    @DisplayName("interpolateString escapes a quote-breakout attempt")
    public void interpolateStringEscapes() {
        FuncParam p = new FuncParam(Types.STRING, "name", "' OR 1=1 --");
        assertEquals("''' OR 1=1 --'", postgres.interpolateString(p));
    }

    @Test
    @DisplayName("addBrackets doubles the closing double-quote (Postgres)")
    public void addBracketsDoublesQuoteChar() {
        assertEquals("\"x\"\" ; drop table t --\"", postgres.addBrackets("x\" ; drop table t --"));
    }

    @Test
    @DisplayName("addBrackets doubles the closing bracket (MSSQL)")
    public void addBracketsDoublesBracketChar() {
        assertEquals("[a]]; drop]", mssql.addBrackets("a]; drop"));
    }

    @Test
    @DisplayName("validateMutationIdentifier rejects a quote/semicolon injection but accepts a plain name")
    public void validateIdentifierRejectsInjection() {
        assertThrows(MutationValidationException.class, () -> postgres.validateMutationIdentifier("public' OR '1'='1"));
        assertThrows(MutationValidationException.class, () -> postgres.validateMutationIdentifier("t; drop table users"));
        assertDoesNotThrow(() -> postgres.validateMutationIdentifier("public"));
        assertDoesNotThrow(() -> postgres.validateMutationIdentifier("my_schema.my_table"));
    }
}
