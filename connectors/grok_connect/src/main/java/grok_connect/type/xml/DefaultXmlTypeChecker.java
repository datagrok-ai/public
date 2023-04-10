package grok_connect.type.xml;

import grok_connect.type.TypeChecker;

public class DefaultXmlTypeChecker implements TypeChecker {
    @Override
    public boolean isSupported(int type, String typeName, int precision, int scale) {
        return (type == java.sql.Types.SQLXML || typeName.equalsIgnoreCase("xml"));
    }
}
