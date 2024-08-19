# Create custom connectors

The Datagrok server uses the [Grok Connect](https://github.com/datagrok-ai/public/blob/master/connectors/README.md) service to connect and retrieve database data. You can extend Grok Connect using a [JVM language](https://www.oracle.com/technical-resources/articles/java/architect-languages.html).

## Adding a new connector

:::info prerequisites

The example on this page uses Java. To follow the instructions, you need the following:

* [Git](https://git-scm.com/) to fetch repository with code
* [Java 8](https://www.java.com/download/ie_manual.jsp) to run and compile code
* [Maven](https://maven.apache.org/download.cgi) to build a project. For languages other than Java, add the appropriate [Maven](https://maven.apache.org/) plugin to 
pom.xml and configure it according to the documentation.
* Code editor of your choice.

:::

To add a new connector:

1. Clone the [Datagrok's public repository](https://github.com/datagrok-ai/public) from GitHub:

   ```bash
   git clone https://github.com/datagrok-ai/public.git
   ```

2. Add a JDBC driver:

   * As a [jar file](https://docs.oracle.com/javase/8/docs/technotes/guides/jar/jarGuide.html) to _public/connectors/grok_connect/src/main/java/grok_connect/lib_ folder.
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

   1. Go to the _providers_ package in the _grok connect_ folder:

      ```bash
      cd public/connectors/grok_connect/src/main/java/grok_connect/providers
      ```

   2. Create a `.java` file and name using the standard naming convention. For example:

      ```
      OrientDbJdbcProvider.java
      ```

   3. In your code editor, open this JAVA file and create a Java class with the same name. Extend it from [JdbcDataProvider](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java), which provides basic functionality for all JDBC providers. 

      For our example, paste this code to your file:

      ```
      public class OrientDbJdbcProvider extends JdbcDataProvider {
          public OrientDbJdbcProvider() {
          }
      }
      ```

      > Note: For simplicity, we omit all imports.

   4. Set the `driverClassName` field using the full driver class name. For our example, within the class constructor:

      ```
      driverClassName = "com.orientechnologies.orient.jdbc.OrientJdbcDriver";
      ```
  
      > Note: To get the driver class name, use the documentation for your chosen JDBC driver.

   5. In the constructor, configure the connection:

      ```
      descriptor = new DataSource();
      descriptor.type = "OrientDb";
      descriptor.description = "Query OrientDb";
      descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
      descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
      descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
      ```

   6. Override the `getConnectionStringImpl` method for your database.

      For our example, add the following code to the class:

      ```
      @Override
      public String getConnectionStringImpl(DataConnection conn) {
          return String.format("jdbc:orient:remote:%s/%s", conn.getServer(), conn.getDb());
      }
      ```

      To complete this step, you need to know how the JDBC connection string is built for your database. Refer to the appropriate documentation.

   7. Register your provider class in [ProviderManager](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/utils/ProviderManager.java) by adding it to the `providersList` in the constructor.

4. Build Grok Connect: 

   Go to the _connectors_ folder with the parent `pom.xml` and run the following command:

   ```bash
   mvn package -DskipTests
   ```

## Testing your connector

:::info prerequisites

A live instance of your chosen database.

:::

Testing options: 

* Using `GrokConnectShell`:

  1. Open `connectors/examples/query.json` and add the necessary details.
  2. Go to the `/public/connectors/grok_connect` folder and run the following command:

      ```bash
      java -cp ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnectShell --q <ABSOLUTE PATH TO query.json>
      ```

* Using Datagrok running locally:

  1. From the `/public/connectors` folder, run:

     ```bash
      java -jar ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnect
     ```
 
  2. Open the platform in a browser and go to **Data** > **Databases**.
  3. In the **Database Manager**, right-click your new connector and select **Add connection...** A parameter dialog opens.
  4. In the parameter dialog, [fill in the connection parameters](databases.md#connecting-to-database) and click **TEST**. Upon successful connection, the database appears in the **Database Manager** under the respective data source. 

  > Note: It may take Datagrok some time to recognize a running Grok Connect. To speed up the process, restart the container.
