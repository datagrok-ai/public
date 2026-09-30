---
title: "Create custom connectors"
description: Extend Grok Connect with a custom Java data connector and JDBC driver to support new database providers.
keywords:
  - custom database connector
  - grok connect
  - jdbc driver java
  - extend grok connect
  - add new data provider
  - java connector development
---

Datagrok server uses 
[Grok Connect](https://github.com/datagrok-ai/grok-connect/blob/master/jdbc/README.md) service to 
query databases. You can extend Grok Connect by developing your own data connectors in Java.

## Adding a new connector

:::info prerequisites

The example on this page uses Java. To follow the instructions, you need the following:

* [Git](https://git-scm.com/) to fetch repository with code
* [Java 8](https://www.java.com/en/download/) to run and compile code
* [Maven](https://maven.apache.org/download.cgi) to build a project. For languages other than Java, add the appropriate [Maven](https://maven.apache.org/) plugin to 
pom.xml and configure it according to the documentation.
* Code editor of your choice.

:::

To add a new connector:

1. Clone the [grok-connect repository](https://github.com/datagrok-ai/grok-connect) from GitHub:

   ```bash
   git clone https://github.com/datagrok-ai/grok-connect.git
   ```

2. Add a JDBC driver:

   * As a [jar file](https://docs.oracle.com/javase/8/docs/technotes/guides/jar/jarGuide.html) to [_grok-connect/jdbc/server/lib_](https://github.com/datagrok-ai/grok-connect/tree/master/jdbc/server/lib) folder.
   * Using `pom.xml` if the driver is available on public repositories. 
  
   For example, let's add the [OrientDB](http://orientdb.org/) connector to Grok Connect. Since it's available on [Maven](https://mvnrepository.com/artifact/com.orientechnologies/orientdb-jdbc), insert the following dependency in the `pom.xml` of the `grok_connect` package:

   ```
   <dependency>
       <groupId>com.orientechnologies</groupId>
       <artifactId>orientdb-jdbc</artifactId>
       <version>3.2.21</version>
   </dependency>
   ```

3. Implement the provider:

   1. Add a new connector class derived from [JdbcDataProvider](https://github.com/datagrok-ai/grok-connect/blob/master/jdbc/server/src/main/java/grok_connect/providers/JdbcDataProvider.java):

      ```
      public class OrientDbJdbcProvider extends JdbcDataProvider {
          public OrientDbJdbcProvider() {
          }
      }
      ```

      > Note: For simplicity, we omit all imports.

   2. Set the `driverClassName` field using the full driver class name. For our example, within the class constructor:

      ```
      driverClassName = "com.orientechnologies.orient.jdbc.OrientJdbcDriver";
      ```
  
      > Note: To get the driver class name, use the documentation for your chosen JDBC driver.

   3. Configure the connection in the constructor:

      ```
      descriptor = new DataSource();
      descriptor.type = "OrientDb";
      descriptor.description = "Query OrientDb";
      descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
      descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
      descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
      ```

   4. Specify how the connection string is built by overriding the `getConnectionStringImpl` method:

      ```
      @Override
      public String getConnectionStringImpl(DataConnection conn) {
          return String.format("jdbc:orient:remote:%s/%s", conn.getServer(), conn.getDb());
      }
      ```

   5. Register your provider class in [ProviderManager](https://github.com/datagrok-ai/grok-connect/blob/master/jdbc/server/src/main/java/grok_connect/utils/ProviderManager.java) by adding it to the `providersList` in the constructor.

4. Build Grok Connect: 

   Go to the [_jdbc_](https://github.com/datagrok-ai/grok-connect/tree/master/jdbc) folder with the parent `pom.xml` and run the following command:

   ```bash
   mvn package -DskipTests
   ```

## Testing your connector

Testing options: 

* Using `GrokConnectShell`:

  1. Open [`jdbc/examples/query.json`](https://github.com/datagrok-ai/grok-connect/blob/master/jdbc/examples/query.json) and add the necessary details.
  2. Go to the [`/grok-connect/jdbc/server`](https://github.com/datagrok-ai/grok-connect/tree/master/jdbc/server) folder and run the following command:

      ```bash
      java -cp ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnectShell --q <ABSOLUTE PATH TO query.json>
      ```

* Using Datagrok running locally:

  1. From the [`/grok-connect/jdbc`](https://github.com/datagrok-ai/grok-connect/tree/master/jdbc) folder, run:

     ```bash
      java -jar ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnect
     ```
 
  2. Open the platform in a browser and go to **Data** > **Databases**.
  3. In the **Database Manager**, right-click your new connector and select **Add connection...** A parameter dialog opens.
  4. In the parameter dialog, [fill in the connection parameters](databases.md#connecting-to-database) and click **TEST**. On successful connection, the database appears in the **Database Manager** under the respective data source.

  > Note: It may take Datagrok some time to recognize a running Grok Connect. To speed up the process, restart the container.
