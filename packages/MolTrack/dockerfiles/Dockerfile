FROM informaticsmatters/rdkit-cartridge-debian:Release_2024_03_6

RUN apt-get update && apt-get install -y python3-venv supervisor curl gosu

RUN python3 -m venv /opt/venv

RUN /opt/venv/bin/pip install --no-cache-dir dg-mol-track==0.3.0

ENV PATH="/opt/venv/bin:$PATH"

ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=postgres
ENV POSTGRES_DB=moltrack

ADD db.sql /docker-entrypoint-initdb.d/01_db.sql
ADD schema.sql /docker-entrypoint-initdb.d/02_schema.sql
ADD setup.sql /docker-entrypoint-initdb.d/03_setup.sql
ADD schema_rdkit.sql /docker-entrypoint-initdb.d/04_schema_rdkit.sql

RUN sed -i 's/\^/\;/g' /docker-entrypoint-initdb.d/04_schema_rdkit.sql

EXPOSE 8000 5432

COPY start.sh /start.sh
RUN sed -i 's/\r$//' start.sh
RUN chmod +x /start.sh

CMD ["/start.sh"]