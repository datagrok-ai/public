FROM datagrok/demo_db_surechembl:2024.1

# Patch base-image CVEs: the 2024-frozen demo parent ships stale Debian
# packages (openssl/gnutls/glibc fixes released since). Refresh them.
RUN apt-get update && apt-get -y upgrade && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV POSTGRES_USER=superuser
ENV POSTGRES_DB=surechembl
ENV DB_NAME=surechembl
ENV DB_PASSWORD=datagrok
ENV DB_USER=datagrok

EXPOSE 5432