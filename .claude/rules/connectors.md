---
paths:
  - connectors/**
---

## GrokConnect Connectors

GrokConnect is a Java/Maven project providing JDBC-based database connectors.

```bash
cd connectors
mvn package -DskipTests          # Build
mvn test                          # Run all tests
mvn test -Dtest=ClassName         # Run specific test class
mvn test -Dtest=ClassName#method  # Run specific test method
java -jar grok_connect/target/grok_connect-*.jar  # Run GrokConnect server locally
```

Source is in `connectors/grok_connect/src/main/java/grok_connect/`.
Tests are in `connectors/grok_connect/src/test/java/grok_connect/`.
