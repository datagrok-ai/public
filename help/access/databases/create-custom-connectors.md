# Create custom connectors

Datagrok server uses a separate service called [Grok Connect](https://github.com/datagrok-ai/public/blob/master/connectors/README.md) to connect and retrieve data from databases. It's easily 
extendable and the only thing you need is a basic knowledge of any popular [JVM language](https://www.oracle.com/technical-resources/articles/java/architect-languages.html) and familiarity with [JDBC](https://download.oracle.com/otndocs/jcp/jdbc-4_2-mrel2-spec/) specification (not required, but nice to have).

## Adding a new connector

For all other languages except Java, you need to add the appropriate [Maven](https://maven.apache.org/) plugin to 
pom.xml and configure it according to the documentation. The example here is focused on using Java programming language.
You need next tools to follow examples:

* [Git](https://git-scm.com/) to fetch repository with code
* [Java 8](https://www.java.com/download/ie_manual.jsp) to run and compile code
* [Maven](https://maven.apache.org/download.cgi) to build project
* Code editor of your choice.

To add a new connector:

* Clone [public](https://github.com/datagrok-ai/public) repository from gitHub:

```bash
git clone https://github.com/datagrok-ai/public.git
```

* Add JDBC driver. You can add driver as a [jar file](https://docs.oracle.com/javase/8/docs/technotes/guides/jar/jarGuide.html) to 
_public/connectors/grok_connect/src/main/java/grok_connect/lib_ folder or by using pom.xml if the driver is available 
on public repositories. For our example, let's say we want to add the [OrientDB](http://orientdb.org/) connector to 
Grok Connect. In our case driver is available on [Maven](https://mvnrepository.com/artifact/com.orientechnologies/orientdb-jdbc), and we can add it to the pom.xml of grok_connect package in the dependency section as follows:

```
<dependency>
    <groupId>com.orientechnologies</groupId>
    <artifactId>orientdb-jdbc</artifactId>
    <version>3.2.21</version>
</dependency>
```

* Navigate to the providers package in grok connect folder:

```bash
cd public/connectors/grok_connect/src/main/java/grok_connect/providers
```

* Create a new file with a .java extension and name it using the standard naming convention used for other providers.
For our example let's call it:

```
 OrientDbJdbcProvider.java
```

* Open created java file in your code editor and create Java class with the same name as the file. 
First, You need to extend it from [JdbcDataProvider](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java) - it provides basic functionality for all JDBC providers.
If you follow example is this tutorial just paste the next block of code to your file:

```
  public class OrientDbJdbcProvider extends JdbcDataProvider {
      public OrientDbJdbcProvider() {
      }
  }
```

>Note: All imports here and in the next blocks of code are omitted for simplicity.

* Override [JdbcDataProvider](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java) field called **driverClassName** and assign it full driver class name.
You can get the name from the appropriate documentation for the JDBC driver of your choice. For our example add the next line to the constructor of the previously created class:

```
  driverClassName = "com.orientechnologies.orient.jdbc.OrientJdbcDriver";
```

Also, paste the next block of code to the constructor. What it does is basically the configuration of the connection form.
  
```
  descriptor = new DataSource();
  descriptor.type = "OrientDb";
  descriptor.description = "Query OrientDb";
  descriptor.connectionTemplate = new ArrayList<>(DbCredentials.dbConnectionTemplate);
  descriptor.connectionTemplate.add(new Property(Property.BOOL_TYPE, DbCredentials.SSL));
  descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate;
```

* Override [JdbcDataProvider](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/providers/JdbcDataProvider.java) method called **getConnectionStringImpl**.
You need to know how JDBC connection string is built for your database. Use appropriate documentation for this.
For our example add the following code to the class:

```
  @Override
    public String getConnectionStringImpl(DataConnection conn) {
        return String.format("jdbc:orient:remote:%s/%s", conn.getServer(), conn.getDb());
    }
```

* Register you provider class in [ProviderManager](https://github.com/datagrok-ai/public/blob/master/connectors/grok_connect/src/main/java/grok_connect/utils/ProviderManager.java). Add it in providersList in the constructor.

* Build Grok Connect. Navigate to _connectors_ folder where parent pom.xml is located and run following command:

```bash
mvn package -DskipTests
```

## Testing added connector

To test your database connector use GrokConnectShell or local build of Datagrok and Grok Connect.

>Note: You need a running instance of the database of your choice to test your newly created connector.

To test with GrokConnectShell:

* Open _connectors/examples/query.json_ and fill it with necessary information
* Navigate to _/public/connectors/grok_connect_ folder and run following command:

```bash
java -cp ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnectShell --q <ABSOLUTE PATH TO query.json>
```

To test with the locally running Datagrok:

* Run from _/public/connectors_ folder:

```bash
 java -jar ./target/<NAME OF GROK CONNECT JAR>.jar grok_connect.GrokConnect
```

>Note: Datagrok needs some time to recognize running Grok Connect or you can reload the container.

* Open platform on the browser and go to **Data** > **Databases** and find connector you have created.
* Right-click the connector, choose **Add connection** and fill the form.
* Test it!
