package grok_connect;

import grok_connect.utils.ConnectionPool;
import grok_connect.utils.Settings;
import grok_connect.utils.SettingsManager;
import org.junit.Test;

import java.sql.Connection;
import java.sql.Statement;
import java.util.Properties;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.CountDownLatch;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;


public class GrokConnectTest {
    ConnectionPool cp = ConnectionPool.getInstance();
    String url = "jdbc:postgresql://localhost:5432/datagrok";
    String driverName = "org.postgresql.Driver";
    String query = "select count(1) from entities";

    @Test
    public void testConnect() {
        //GrokConnect.connect("test");
    }

//    @Test
//    public void connectionPoolTest() {
//        SettingsManager.getInstance().initSettingsWithDefaults();
//
//        Properties prop1 = new Properties();
//        prop1.setProperty("user", "datagrok");
//        prop1.setProperty("password", "datagrok");
//        try {
//            cp.getConnection(url, prop1, driverName).createStatement(); // create connection
//        } catch (Exception throwables) {
//            fail();
//        }
//
//        Properties prop2 = new Properties(); // different user
//        prop2.setProperty("user", "test_user");
//        prop2.setProperty("password", "test_user");
//        try {
//            cp.getConnection(url, prop2, driverName).createStatement(); // create connection
//        } catch (Exception throwables) {
//            fail();
//        }
//
//        try {
//            Statement statement = cp.getConnection(url, prop1, driverName).createStatement(); // use connection
//            statement.execute(query);
//        } catch (Exception throwables) {
//            fail();
//        }
//
//        try {
//            Statement statement = cp.getConnection(url, prop2, driverName).createStatement(); // use connection
//            statement.execute(query);
//            fail();
//        } catch (Exception throwables) {
//            assertEquals("ERROR: permission denied for table entities", throwables.getMessage());
//        }
//
//        assertEquals(2, cp.connectionPool.size());
//    }
//
//
//    @Test
//    public void connectionPoolTimeoutTest() throws InterruptedException {
//        SettingsManager.getInstance().initSettingsWithDefaults();
//        SettingsManager.getInstance().settings.connectionPoolIdleTimeout = 10000;
//
//        Properties prop = new Properties();
//        prop.setProperty("user", "datagrok");
//        prop.setProperty("password", "datagrok");
//
//        try {
//            Connection conn = cp.getConnection(url, prop, driverName);
//            Statement statement = conn.createStatement();
//            statement.execute(query);
//            conn.close();
//        } catch (Exception throwables) {
//            fail();
//        }
//
//        Thread.sleep(SettingsManager.getInstance().settings.connectionPoolIdleTimeout + 30000); //plus maximum variation
//        String key = url + prop + driverName;
//        assertEquals(0, cp.connectionPool.get(key).poolProxy.getTotalConnections());
//    }
}