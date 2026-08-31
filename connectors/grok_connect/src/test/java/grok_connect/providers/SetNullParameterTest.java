package grok_connect.providers;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;
import serialization.Types;

import java.sql.PreparedStatement;
import java.sql.SQLException;

import static org.mockito.Mockito.verify;

public class SetNullParameterTest {
    private final PostgresDataProvider provider = new PostgresDataProvider();

    @Test
    @DisplayName("null bigint binds as SQL BIGINT, not VARCHAR")
    public void bigIntNullBindsAsBigInt() throws SQLException {
        PreparedStatement st = Mockito.mock(PreparedStatement.class);
        provider.setNullParameter(st, 1, Types.BIG_INT);
        verify(st).setNull(1, java.sql.Types.BIGINT);
    }

    @Test
    @DisplayName("null bool binds as SQL BOOLEAN, not VARCHAR")
    public void boolNullBindsAsBoolean() throws SQLException {
        PreparedStatement st = Mockito.mock(PreparedStatement.class);
        provider.setNullParameter(st, 1, Types.BOOL);
        verify(st).setNull(1, java.sql.Types.BOOLEAN);
    }

    @Test
    @DisplayName("null int/double bind as SQL NUMERIC")
    public void numericNullBindsAsNumeric() throws SQLException {
        PreparedStatement st = Mockito.mock(PreparedStatement.class);
        provider.setNullParameter(st, 2, Types.INT);
        verify(st).setNull(2, java.sql.Types.NUMERIC);
        provider.setNullParameter(st, 3, Types.FLOAT);
        verify(st).setNull(3, java.sql.Types.NUMERIC);
    }

    @Test
    @DisplayName("null string falls back to SQL VARCHAR")
    public void stringNullBindsAsVarchar() throws SQLException {
        PreparedStatement st = Mockito.mock(PreparedStatement.class);
        provider.setNullParameter(st, 1, Types.STRING);
        verify(st).setNull(1, java.sql.Types.VARCHAR);
    }
}
