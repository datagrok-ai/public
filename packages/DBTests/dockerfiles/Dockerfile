FROM postgres
ENV POSTGRES_PASSWORD datagrok
ENV POSTGRES_USER datagrok
ENV POSTGRES_DB world
COPY world.sql /docker-entrypoint-initdb.d/
EXPOSE 5432

HEALTHCHECK --interval=5s --timeout=3s --retries=10 --start-period=10s \
  CMD pg_isready -U "$POSTGRES_USER" -d "$POSTGRES_DB" || exit 1