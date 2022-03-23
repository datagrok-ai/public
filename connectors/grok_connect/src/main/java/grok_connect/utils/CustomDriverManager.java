/*
 * Copyright (c) 1996, 2013, Oracle and/or its affiliates. All rights reserved.
 * ORACLE PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

package grok_connect.utils;
import java.sql.*;
import java.util.*;
import java.security.AccessController;
import java.security.PrivilegedAction;

//import sun.reflect.CallerSensitive;


/**
 * <P>The basic service for managing a set of JDBC drivers.<br>
 * <B>NOTE:</B> The {@link javax.sql.DataSource} interface, new in the
 * JDBC 2.0 API, provides another way to connect to a data source.
 * The use of a <code>DataSource</code> object is the preferred means of
 * connecting to a data source.
 *
 * <P>As part of its initialization, the <code>DriverManager</code> class will
 * attempt to load the driver classes referenced in the "jdbc.drivers"
 * system property. This allows a user to customize the JDBC Drivers
 * used by their applications. For example in your
 * ~/.hotjava/properties file you might specify:
 * <pre>
 * <CODE>jdbc.drivers=foo.bah.Driver:wombat.sql.Driver:bad.taste.ourDriver</CODE>
 * </pre>
 *<P> The <code>DriverManager</code> methods <code>getConnection</code> and
 * <code>getDrivers</code> have been enhanced to support the Java Standard Edition
 * <a href="../../../technotes/guides/jar/jar.html#Service%20Provider">Service Provider</a> mechanism. JDBC 4.0 Drivers must
 * include the file <code>META-INF/services/java.sql.Driver</code>. This file contains the name of the JDBC drivers
 * implementation of <code>java.sql.Driver</code>.  For example, to load the <code>my.sql.Driver</code> class,
 * the <code>META-INF/services/java.sql.Driver</code> file would contain the entry:
 * <pre>
 * <code>my.sql.Driver</code>
 * </pre>
 *
 * <P>Applications no longer need to explicitly load JDBC drivers using <code>Class.forName()</code>. Existing programs
 * which currently load JDBC drivers using <code>Class.forName()</code> will continue to work without
 * modification.
 *
 * <P>When the method <code>getConnection</code> is called,
 * the <code>DriverManager</code> will attempt to
 * locate a suitable driver from amongst those loaded at
 * initialization and those loaded explicitly using the same classloader
 * as the current applet or application.
 *
 * <P>
 * Starting with the Java 2 SDK, Standard Edition, version 1.3, a
 * logging stream can be set only if the proper
 * permission has been granted.  Normally this will be done with
 * the tool PolicyTool, which can be used to grant <code>permission
 * java.sql.SQLPermission "setLog"</code>.
 * @see Driver
 * @see Connection
 */
public class CustomDriverManager {


    // List of registered JDBC drivers
    private static volatile int loginTimeout = 0;
    private static volatile java.io.PrintWriter logWriter = null;
    private static volatile java.io.PrintStream logStream = null;
    // Used in println() to synchronize logWriter
    private final static  Object logSync = new Object();

    /* Prevent the DriverManager class from being instantiated. */
    private CustomDriverManager(){}


    /**
     * Load the initial JDBC drivers by checking the System property
     * jdbc.properties and then use the {@code ServiceLoader} mechanism
     */
    static {
        loadInitialDrivers();
        println("JDBC DriverManager initialized");
    }

    /**
     * The <code>SQLPermission</code> constant that allows the
     * setting of the logging stream.
     * @since 1.3
     */
    final static SQLPermission SET_LOG_PERMISSION =
            new SQLPermission("setLog");

    /**
     * The {@code SQLPermission} constant that allows the
     * un-register a registered JDBC driver.
     * @since 1.8
     */
    final static SQLPermission DEREGISTER_DRIVER_PERMISSION =
            new SQLPermission("deregisterDriver");

    //--------------------------JDBC 2.0-----------------------------

    //---------------------------------------------------------------

    /**
     * Sets the logging/tracing PrintStream that is used
     * by the <code>DriverManager</code>
     * and all drivers.
     *<P>
     * In the Java 2 SDK, Standard Edition, version 1.3 release, this method checks
     * to see that there is an <code>SQLPermission</code> object before setting
     * the logging stream.  If a <code>SecurityManager</code> exists and its
     * <code>checkPermission</code> method denies setting the log writer, this
     * method throws a <code>java.lang.SecurityException</code>.
     *
     * @param out the new logging/tracing PrintStream; to disable, set to <code>null</code>
     * @deprecated Use {@code setLogWriter}
     * @throws SecurityException if a security manager exists and its
     *    <code>checkPermission</code> method denies setting the log stream
     *
     * @see SecurityManager#checkPermission
     * @see #getLogStream
     */
    @Deprecated
    public static void setLogStream(java.io.PrintStream out) {

        SecurityManager sec = System.getSecurityManager();
        if (sec != null) {
            sec.checkPermission(SET_LOG_PERMISSION);
        }

        logStream = out;
        if ( out != null )
            logWriter = new java.io.PrintWriter(out);
        else
            logWriter = null;
    }

    /**
     * Retrieves the logging/tracing PrintStream that is used by the <code>DriverManager</code>
     * and all drivers.
     *
     * @return the logging/tracing PrintStream; if disabled, is <code>null</code>
     * @deprecated  Use {@code getLogWriter}
     * @see #setLogStream
     */
    @Deprecated
    public static java.io.PrintStream getLogStream() {
        return logStream;
    }

    /**
     * Prints a message to the current JDBC log stream.
     *
     * @param message a log or tracing message
     */
    public static void println(String message) {
        synchronized (logSync) {
            if (logWriter != null) {
                logWriter.println(message);

                // automatic flushing is never enabled, so we must do it ourselves
                logWriter.flush();
            }
        }
    }

    //------------------------------------------------------------------------

    // Indicates whether the class object that would be created if the code calling
    // DriverManager is accessible.
    private static boolean isDriverAllowed(Driver driver, Class<?> caller) {
        ClassLoader callerCL = caller != null ? caller.getClassLoader() : null;
        return isDriverAllowed(driver, callerCL);
    }

    private static boolean isDriverAllowed(Driver driver, ClassLoader classLoader) {
        boolean result = false;
        if(driver != null) {
            Class<?> aClass = null;
            try {
                aClass =  Class.forName(driver.getClass().getName(), true, classLoader);
            } catch (Exception ex) {
                result = false;
            }

            result = ( aClass == driver.getClass() ) ? true : false;
        }

        return result;
    }

    private static void loadInitialDrivers() {
        String drivers;
        try {
            drivers = AccessController.doPrivileged(new PrivilegedAction<String>() {
                public String run() {
                    return System.getProperty("jdbc.drivers");
                }
            });
        } catch (Exception ex) {
            drivers = null;
        }
        // If the driver is packaged as a Service Provider, load it.
        // Get all the drivers through the classloader
        // exposed as a java.sql.Driver.class service.
        // ServiceLoader.load() replaces the sun.misc.Providers()

        AccessController.doPrivileged(new PrivilegedAction<Void>() {
            public Void run() {

                ServiceLoader<Driver> loadedDrivers = ServiceLoader.load(Driver.class);
                Iterator<Driver> driversIterator = loadedDrivers.iterator();

                /* Load these drivers, so that they can be instantiated.
                 * It may be the case that the driver class may not be there
                 * i.e. there may be a packaged driver with the service class
                 * as implementation of java.sql.Driver but the actual class
                 * may be missing. In that case a java.util.ServiceConfigurationError
                 * will be thrown at runtime by the VM trying to locate
                 * and load the service.
                 *
                 * Adding a try catch block to catch those runtime errors
                 * if driver not available in classpath but it's
                 * packaged as service and that service is there in classpath.
                 */
                try{
                    while(driversIterator.hasNext()) {
                        driversIterator.next();
                    }
                } catch(Throwable t) {
                    // Do nothing
                }
                return null;
            }
        });

        println("DriverManager.initialize: jdbc.drivers = " + drivers);

        if (drivers == null || drivers.equals("")) {
            return;
        }
        String[] driversList = drivers.split(":");
        println("number of Drivers:" + driversList.length);
        for (String aDriver : driversList) {
            try {
                println("DriverManager.Initialize: loading " + aDriver);
                Class.forName(aDriver, true,
                        ClassLoader.getSystemClassLoader());
            } catch (Exception ex) {
                println("DriverManager.Initialize: load failed: " + ex);
            }
        }
    }

    public static void printDriverNames() {
        for(Driver aDriver : Collections.list(DriverManager.getDrivers())) {
            System.out.println(aDriver.getClass().getName());
        }
    }

    //  Worker method called by the public getConnection() methods.
    public static Connection getConnection(
            String url, java.util.Properties info, String driverClassName) throws SQLException {
//        /*
//         * When callerCl is null, we should check the application's
//         * (which is invoking this class indirectly)
//         * classloader, so that the JDBC driver class outside rt.jar
//         * can be loaded from here.
//         */
//        ClassLoader callerCL = caller != null ? caller.getClassLoader() : null;
//        synchronized(DriverManager.class) {
//            // synchronize loading of the correct classloader.
//            if (callerCL == null) {
//                callerCL = Thread.currentThread().getContextClassLoader();
//            }
//        }

        if(url == null || url.equals("")) {
            throw new SQLException("The url cannot be null", "08001");
        }

        Connection conn = ConnectionPool.getInstance().getNativeConnection(url, info, driverClassName);
        if (conn != null && !conn.isClosed())
            return conn;

        println("DriverManager.getConnection(\"" + url + "\")");

        // Walk through the loaded registeredDrivers attempting to make a connection.
        // Remember the first exception that gets raised so we can reraise it.
        SQLException reason = null;

        for(Driver aDriver : Collections.list(DriverManager.getDrivers())) {
            // If the caller does not have permission to load the driver then
            // skip it.
            if(aDriver.getClass().getName().equals(driverClassName)) {
                try {
                    println("    trying " + aDriver.getClass().getName());
                    Connection con = aDriver.connect(url, info);
                    if (con != null) {
                        // Success!
                        println("getConnection returning " + aDriver.getClass().getName());
                        ConnectionPool.getInstance().putNativeConnection(con, url, info, driverClassName);
                        return (con);
                    }
                } catch (SQLException ex) {
                    if (reason == null) {
                        reason = ex;
                    }
                }

            } else {
                println("    skipping: " + aDriver.getClass().getName());
            }

        }

        // if we got here nobody could connect.
        if (reason != null)    {
            println("getConnection failed: " + reason);
            throw reason;
        }

        println("getConnection: no suitable driver found for "+ url);
        throw new SQLException("No suitable driver found for "+ url, "08001");
    }


}
