package grok_connect;

import grok_connect.utils.ConnectionPool;
import org.junit.Test;

import java.sql.SQLException;
import java.sql.Statement;
import java.util.Properties;

import static org.junit.Assert.*;


public class GrokConnectTest {
    @Test
    public void testConnect() {
        //GrokConnect.connect("test");
    }

    @Test
    public void connectionPoolTest() {
        ConnectionPool cp = ConnectionPool.getInstance();
        String url = "jdbc:postgresql://localhost:5432/datagrok";
        String driverName = "org.postgresql.Driver";
        String query = "select count(1) from entities";

        Properties prop1 = new Properties();
        prop1.setProperty("user", "datagrok");
        prop1.setProperty("password", "datagrok");
        try {
            cp.getConnection(url, prop1, driverName).createStatement(); // create connection
        } catch (SQLException throwables) {
            fail();
        }

        Properties prop2 = new Properties(); // different user
        prop2.setProperty("user", "test_user");
        prop2.setProperty("password", "test_user");
        try {
            cp.getConnection(url, prop2, driverName).createStatement(); // create connection
        } catch (SQLException throwables) {
            fail();
        }

        try {
            Statement statement = cp.getConnection(url, prop1, driverName).createStatement(); // use connection
            statement.execute(query);
        } catch (SQLException throwables) {
            fail();
        }

        try {
            Statement statement = cp.getConnection(url, prop2, driverName).createStatement(); // use connection
            statement.execute(query);
            fail();
        } catch (SQLException throwables) {
            assertEquals("ERROR: permission denied for table entities", throwables.getMessage());
        }

        assertEquals(2, cp.connectionPool.size());
    }
}