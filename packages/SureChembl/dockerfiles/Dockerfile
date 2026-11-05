FROM datagrok/demo_db_surechembl:2024.1 AS src

# The 2024-frozen demo parent is Ubuntu 22.04 with hundreds of stale OS packages. Only the
# self-contained conda env with PostgreSQL 16 + RDKit cartridge (/opt/rdkit), the prebuilt database
# and the entrypoint scripts are carried over onto Wolfi (glibc, near-zero CVEs). The unused
# /opt/py27 data-prep env (Python 2.7) is left behind.
FROM cgr.dev/chainguard/wolfi-base

RUN apk upgrade --no-cache && apk add --no-cache bash glibc-locale-en && \
    addgroup -g 2001 grok && \
    adduser -D -u 1001 -G grok -h /home/grok -s /bin/bash grok

COPY --from=src /opt/rdkit /opt/rdkit
COPY --from=src /usr/local/bin/ /usr/local/bin/
COPY --from=src /docker-entrypoint-initdb.d /docker-entrypoint-initdb.d
COPY --from=src --chown=1001:2001 /home/grok /home/grok

ENV PATH=/opt/rdkit/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
ENV PGDATA=/home/grok/postgresql/data
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
ENV DB_SCHEMA=public,rdk
ENV POSTGRES_USER=superuser
ENV POSTGRES_DB=surechembl
ENV DB_NAME=surechembl
ENV DB_PASSWORD=datagrok
ENV DB_USER=datagrok

USER grok
# Proves the conda PostgreSQL and the RDKit cartridge run on this base: the healthcheck's query on a scratch copy.
RUN cp -a "$PGDATA" /tmp/pgcheck && \
    pg_ctl -D /tmp/pgcheck -o "-c listen_addresses='' -k /tmp" -w start && \
    test "$(psql -h /tmp -U "$POSTGRES_USER" -d "$POSTGRES_DB" -tAc "select count(*) from rdk.mols where m@>'CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1'")" = 31 && \
    pg_ctl -D /tmp/pgcheck -m fast -w stop && \
    rm -rf /tmp/pgcheck

VOLUME /home/grok/postgresql/data
STOPSIGNAL SIGINT
HEALTHCHECK CMD ["docker-healthcheck.sh"]
EXPOSE 5432
ENTRYPOINT ["entrypoint.sh"]
CMD ["postgres"]
