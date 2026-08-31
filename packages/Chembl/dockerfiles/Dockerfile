# The ChEMBL release this package ships is the base image tag, and the package
# version is that release plus `.0.0` — bump both together.
FROM datagrok/demo_db_chembl:37

# The base only sets the superuser; the demo-db entrypoint needs these to create
# the read-only role the package connections log in as.
ENV DB_NAME=chembl
ENV DB_USER=datagrok
ENV DB_PASSWORD=datagrok

EXPOSE 5432
