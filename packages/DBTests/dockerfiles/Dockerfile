FROM postgres:17-alpine

# su-exec replaces gosu: the base image's gosu is a Go binary whose embedded
# toolchain carries CRITICAL/HIGH CVEs on every scan; su-exec has the same CLI.
# Deleted, not overwritten: scanners keep reporting a lower-layer binary replaced in place.
# The entrypoint finds gosu on PATH, where /usr/local/sbin comes first.
RUN apk upgrade --no-cache && apk add --no-cache su-exec && \
    rm /usr/local/bin/gosu && mkdir -p /usr/local/sbin && ln -s /sbin/su-exec /usr/local/sbin/gosu && \
    gosu postgres true

ENV POSTGRES_PASSWORD datagrok
ENV POSTGRES_USER datagrok
ENV POSTGRES_DB world
# musl has no locale collations: use ICU en-US to keep the en_US.utf8 sort order of the Debian image
ENV POSTGRES_INITDB_ARGS="--locale-provider=icu --icu-locale=en-US"
COPY world.sql /docker-entrypoint-initdb.d/
EXPOSE 5432

HEALTHCHECK --interval=5s --timeout=3s --retries=10 --start-period=10s \
  CMD pg_isready -U "$POSTGRES_USER" -d "$POSTGRES_DB" || exit 1
