FROM postgres
ENV POSTGRES_PASSWORD datagrok
ENV POSTGRES_USER datagrok
ENV POSTGRES_DB world
COPY world.sql /docker-entrypoint-initdb.d/
EXPOSE 5432
