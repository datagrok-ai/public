package grok_connect.providers

import grok_connect.connectors_info.DataConnection
import grok_connect.connectors_info.DataSource
import grok_connect.connectors_info.DbCredentials
import grok_connect.resultset.ResultSetManager
import grok_connect.table_query.AggrFunctionInfo
import grok_connect.table_query.Stats
import grok_connect.utils.Prop
import grok_connect.utils.Property
import grok_connect.utils.ProviderManager
import serialization.Types
import utilities.extetions.orIfNull
import utilities.extetions.substringOrEmpty
import utilities.extetions.takeIfTrue

private const val SPACE = " "
private const val FIRST_INDEX = 0
private const val FIELD_SCHEMA_NAME = "schema_name"
private const val FIELD_ALIAS_SCHEMA_NAME = "table_schema"

class SapHanaDataProvider(resultSetManager:ResultSetManager, providerManager: ProviderManager) : JdbcDataProvider(resultSetManager, providerManager) {

    init {
        driverClassName = "com.sap.db.jdbc.Driver"
        descriptor = DataSource()
        descriptor.type = "SAP HANA"
        descriptor.description = "Query SAP HANA database"
        descriptor.connectionTemplate = getConnectionTemplate()
        descriptor.credentialsTemplate = DbCredentials.dbCredentialsTemplate
        descriptor.canBrowseSchema = true
        descriptor.nameBrackets = "\""
        descriptor.aggregations.add(AggrFunctionInfo(Stats.STDEV, "stddev(#)", Types.dataFrameNumericTypes))

    }

    override fun getConnectionStringImpl(conn: DataConnection): String {
        val port = conn.port?.let { ":${conn.port}" }.orEmpty()

        return buildString{
            append("jdbc:sap://")
            append(conn.server)
            append(port)
            append("/?databaseName=${conn.db}")
        }

    }

    override fun getSchemasSql(db: String?): String {
        val sysSchemasFilter = getSysSchemasFilter()

        return buildString{
            append("SELECT $FIELD_SCHEMA_NAME AS $FIELD_ALIAS_SCHEMA_NAME FROM")
            append(SPACE)
            append("(")
            append("SELECT DISTINCT $FIELD_SCHEMA_NAME FROM tables")
            append(SPACE)
            append("UNION")
            append(SPACE)
            append("SELECT DISTINCT $FIELD_SCHEMA_NAME FROM views")
            append(")")
            append(SPACE)
            append("WHERE $sysSchemasFilter")
            append(SPACE)
            append("ORDER BY $FIELD_ALIAS_SCHEMA_NAME")
        }

    }

    override fun getSchemaSql(db: String?, schema: String?, table: String?): String {
        val sysSchemasFilter = getSysSchemasFilter()
        var whereClause = "WHERE $sysSchemasFilter"
        table?.also { whereClause = "$whereClause AND table_name = '$table'" }
        schema?.also { whereClause = "$whereClause AND $FIELD_SCHEMA_NAME = '$schema'" }

        return buildString{
            append("SELECT $FIELD_SCHEMA_NAME AS $FIELD_ALIAS_SCHEMA_NAME, table_name, column_name, data_type_name AS data_type")
            append(SPACE)
            append("FROM table_columns")
            append(SPACE)
            append(whereClause)
            append(SPACE)
            append("ORDER BY table_name")
        }

    }

    override fun limitToSql(query: String, limit: Int): String {
        return "SELECT TOP $limit * FROM ($query)"

    }

    override fun addBrackets(name: String): String {
        val brackets = descriptor.nameBrackets
        val nameStartsWith = name.startsWith(brackets.substringOrEmpty(0, 1))
        val startBracket = brackets.substringOrEmpty(0, 1)
        val endBracket = brackets.substringOrEmpty(brackets.length - 1, brackets.length)

        return nameStartsWith
            .takeIfTrue()
            ?.let { name }
            .orIfNull { "$startBracket$name$endBracket" }

    }

    private fun getConnectionTemplate(): MutableList<Property> {
        return mutableListOf(
            Property(Property.STRING_TYPE, DbCredentials.SERVER, DbCredentials.DB_DESCRIPTION),
            Property(Property.INT_TYPE, DbCredentials.PORT, Prop()),
            Property(Property.STRING_TYPE, DbCredentials.DB, DbCredentials.DB_DESCRIPTION)
        )

    }

    private fun getSysSchemasFilter(): String {
        return SysSchemasFilter
            .values()
            .mapIndexed { index, sysSchemasFilter ->
                val schema = sysSchemasFilter.value
                val condition = "$FIELD_SCHEMA_NAME LIKE '$schema'"
                index.takeIf { it == FIRST_INDEX }
                    ?.let { "NOT $condition" }
                    .orIfNull { "AND NOT $condition"  }
            }.joinToString(separator = SPACE)

    }

}

private enum class SysSchemasFilter(val value: String) {
    FDT("%FDT%"),
    SAP("%SAP%"),
    SYS("%SYS%"),
    UIS("UIS"),
    HANA_XS_BASE("HANA_XS_BASE")
}
