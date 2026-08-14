package grok_connect.utils;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

/**
 * Classification matrix for the §6.2 read-only query policy (connector-writes WO-B13):
 * the classifier is a read-WHITELIST (SELECT / WITH / SHOW / EXPLAIN / DESCRIBE / VALUES) —
 * anything else, including empty or unparseable statements, classifies as a write.
 */
class StatementClassifierTest {

    @DisplayName("Read shapes: whitelisted first keyword after whitespace/comments/parens")
    @ParameterizedTest
    @ValueSource(strings = {
            "SELECT 1",
            "  select * from t",
            "\tSELECT\n1",
            "-- leading comment\nSELECT 1",
            "-- a\n-- b\nselect 1",
            "/* block */ SELECT 1",
            "/* multi\nline */\nSELECT 1",
            "(SELECT 1) UNION (SELECT 2)",
            ";;SELECT 1",
            "WITH t AS (SELECT 1) SELECT * FROM t",
            "SHOW server_version",
            "EXPLAIN SELECT 1",
            "DESCRIBE t",
            "VALUES (1)"
    })
    public void readShapes(String statement) {
        Assertions.assertTrue(StatementClassifier.isRead(statement, "--"),
                "expected READ: " + statement);
    }

    @DisplayName("Write shapes: any non-whitelisted first keyword classifies as a write")
    @ParameterizedTest
    @ValueSource(strings = {
            "INSERT INTO t VALUES (1)",
            "update t set a = 1",
            "DELETE FROM t",
            "MERGE INTO t USING s ON (t.id = s.id) WHEN MATCHED THEN UPDATE SET a = 1",
            "CREATE TABLE t (id int)",
            "DROP TABLE t",
            "ALTER TABLE t ADD COLUMN c int",
            "TRUNCATE t",
            "GRANT ALL ON t TO u",
            "CALL p()",
            "COPY t FROM stdin",
            "VACUUM",
            "SET role admin",
            "BEGIN",
            "DO $$ BEGIN NULL; END $$",
            "-- comment then write\nINSERT INTO t VALUES (1)",
            "/* block */ DROP TABLE t"
    })
    public void writeShapes(String statement) {
        Assertions.assertFalse(StatementClassifier.isRead(statement, "--"),
                "expected WRITE: " + statement);
    }

    @DisplayName("Degenerate shapes classify as writes (conservative refuse)")
    @Test
    public void degenerateShapes() {
        Assertions.assertFalse(StatementClassifier.isRead(null, "--"));
        Assertions.assertFalse(StatementClassifier.isRead("", "--"));
        Assertions.assertFalse(StatementClassifier.isRead("   ", "--"));
        Assertions.assertFalse(StatementClassifier.isRead("-- only a comment", "--"));
        Assertions.assertFalse(StatementClassifier.isRead("/* only a block */", "--"));
        Assertions.assertFalse(StatementClassifier.isRead("/* unterminated SELECT 1", "--"));
        Assertions.assertFalse(StatementClassifier.isRead("1 + 1", "--"));
    }

    @DisplayName("Provider comment marker is honored (# line comments)")
    @Test
    public void providerCommentMarker() {
        Assertions.assertTrue(StatementClassifier.isRead("# note\nSELECT 1", "#"));
        // with a '--' marker, '#' is not a comment: first token is not a letter keyword -> write
        Assertions.assertFalse(StatementClassifier.isRead("# note\nSELECT 1", "--"));
    }

    @DisplayName("firstKeyword extraction")
    @Test
    public void firstKeywordExtraction() {
        Assertions.assertEquals("SELECT", StatementClassifier.firstKeyword("-- c\n  (select 1)", "--"));
        Assertions.assertEquals("INSERT", StatementClassifier.firstKeyword("/* c */ insert into t", "--"));
        Assertions.assertNull(StatementClassifier.firstKeyword("-- only", "--"));
        Assertions.assertNull(StatementClassifier.firstKeyword(null, "--"));
    }
}
