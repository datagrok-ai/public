FROM datagrok/demo_db_surechembl:2024.1

# Patch base-image CVEs: the 2024-frozen demo parent ships stale Debian
# packages (openssl/gnutls/glibc fixes released since). Refresh them.
# The base runs as the unprivileged grok user; apt needs root, then revert.
# Also drop the unused /opt/py27 data-prep env (Python 2.7): its stale
# pip/setuptools/wheel/certifi have no py2.7 upgrades and are not used at
# runtime (the container only serves postgres). /opt/rdkit is kept intact.
USER root
RUN apt-get update && apt-get -y upgrade && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    rm -rf /opt/py27
USER grok

ENV POSTGRES_USER=superuser
ENV POSTGRES_DB=surechembl
ENV DB_NAME=surechembl
ENV DB_PASSWORD=datagrok
ENV DB_USER=datagrok

EXPOSE 5432