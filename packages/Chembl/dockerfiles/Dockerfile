# The ChEMBL release this package ships is the base image tag, and the package
# version is that release plus `.0.0` — bump both together.
FROM datagrok/demo_db_chembl:37

# The base only sets the superuser; the demo-db entrypoint needs these to create
# the read-only role the package connections log in as.
ENV DB_NAME=chembl
ENV DB_USER=datagrok
ENV DB_PASSWORD=datagrok

# The base's check connects to `hostname -i`, which lists every address of a container on
# more than one network (the spawner attaches two), so it never passes.
HEALTHCHECK CMD pg_isready -h 127.0.0.1 -p 5432

EXPOSE 5432
