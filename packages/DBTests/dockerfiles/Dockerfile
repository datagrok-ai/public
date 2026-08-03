# su-exec replaces gosu: the base image's gosu is a Go binary whose embedded
# toolchain (go 1.24.x) carries CRITICAL/HIGH CVEs on every scan; su-exec is
# ~100 lines of C with the same CLI (su-exec user cmd == gosu user cmd).
FROM postgres:17-bookworm AS suexec
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev wget ca-certificates && \
    wget -qO /tmp/su-exec.c https://raw.githubusercontent.com/ncopa/su-exec/master/su-exec.c && \
    gcc -Wall -Werror -g -o /usr/local/bin/su-exec /tmp/su-exec.c

FROM postgres:17-bookworm

# Patch base-image CVEs: the postgres base does not apt-upgrade its Debian
# packages, so gnutls28 / perl / glibc etc. ship stale. Refresh them.
RUN apt-get update && apt-get -y upgrade && rm -rf /var/lib/apt/lists/*

COPY --from=suexec /usr/local/bin/su-exec /usr/local/bin/gosu

ENV POSTGRES_PASSWORD datagrok
ENV POSTGRES_USER datagrok
ENV POSTGRES_DB world
COPY world.sql /docker-entrypoint-initdb.d/
EXPOSE 5432

HEALTHCHECK --interval=5s --timeout=3s --retries=10 --start-period=10s \
  CMD pg_isready -U "$POSTGRES_USER" -d "$POSTGRES_DB" || exit 1
