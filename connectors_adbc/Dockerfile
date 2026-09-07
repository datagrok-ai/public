# Multi-stage build: ADBC drivers + grok_connect_adbc builder + runtime.

# ---------------------------------------------------------------------------
# Stage 1: build the ClickHouse ADBC driver (libadbc_clickhouse.so) from source.
# ---------------------------------------------------------------------------
FROM rust:1.91-bookworm AS clickhouse-driver

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    protobuf-compiler \
    pkg-config \
    libssl-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Pin the driver revision. Public repo, anonymous HTTPS clone.
ARG ADBC_CLICKHOUSE_REPO=https://github.com/ClickHouse/adbc_clickhouse
ARG ADBC_CLICKHOUSE_REV=d1cc2ead3e1b05ace196d524bff8a91b51345d61

RUN git clone "${ADBC_CLICKHOUSE_REPO}" /src \
    && git -C /src checkout "${ADBC_CLICKHOUSE_REV}"

WORKDIR /src

# `--features ffi` exports the C-ABI entrypoint `AdbcClickhouseInit`; the
# cdylib crate-type produces target/release/libadbc_clickhouse.so.
RUN cargo build --release --features ffi

# ---------------------------------------------------------------------------
# Stage 2: build the BigQuery and Snowflake ADBC drivers from source.
#
# Both are Go builds. Upstream ships them as a PyPI wheel (BigQuery) and an Apache Arrow
# apt package (Snowflake), but those binaries are stamped with whatever Go toolchain and
# dependency set the release was cut with, and ADBC 24 dropped both drivers — so there is
# no newer upstream build to move to. Building here lets us pick the toolchain and raise
# the vulnerable dependency floors; the ldflags mirror the upstream pkg/Makefile.
# ---------------------------------------------------------------------------
FROM golang:1.26-bookworm AS adbc-go-drivers

ARG ADBC_GO_REF=go/adbc/v1.11.0
ARG ADBC_DRIVER_VERSION=v1.11.0

RUN git clone --depth 1 --branch "${ADBC_GO_REF}" https://github.com/apache/arrow-adbc /adbc

WORKDIR /adbc/go/adbc

RUN go get \
        golang.org/x/crypto@v0.54.0 \
        golang.org/x/net@v0.57.0 \
        golang.org/x/text@v0.40.0 \
        github.com/apache/thrift@v0.24.0 \
        github.com/klauspost/compress@v1.19.1 \
        google.golang.org/grpc@v1.83.0 \
        go.opentelemetry.io/otel@v1.45.0 \
        go.opentelemetry.io/otel/sdk@v1.45.0 \
        go.opentelemetry.io/otel/exporters/otlp/otlptrace@v1.45.0 \
        go.opentelemetry.io/otel/exporters/otlp/otlptrace/otlptracegrpc@v1.45.0 \
        go.opentelemetry.io/otel/exporters/otlp/otlptrace/otlptracehttp@v1.45.0 \
        go.opentelemetry.io/otel/exporters/stdout/stdouttrace@v1.45.0 \
        github.com/aws/aws-sdk-go-v2/aws/protocol/eventstream@v1.7.16 \
        github.com/aws/aws-sdk-go-v2/service/s3@v1.106.4 \
    && go mod tidy

RUN mkdir -p /drivers && cd pkg \
    && for d in bigquery snowflake; do \
         go build -buildvcs=false -tags driverlib -buildmode=c-shared \
           -ldflags "-s -w -X github.com/apache/arrow-adbc/go/adbc/driver/internal/driverbase.infoDriverVersion=${ADBC_DRIVER_VERSION}" \
           -o "/drivers/libadbc_driver_${d}.so" "./${d}"; \
       done \
    && rm -f /drivers/*.h

# ---------------------------------------------------------------------------
# Stage 3: build grok_connect_adbc.
# ---------------------------------------------------------------------------
FROM rust:1.91-bookworm AS builder

# Build dependencies:
# - cmake, protobuf-compiler: Arrow / gRPC crates
# - git:                      git dependencies (databricks-adbc)
RUN apt-get update && apt-get install -y \
    cmake \
    protobuf-compiler \
    curl \
    git \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# The databricks-adbc git dep has a C# submodule whose .gitmodules uses an SSH URL.
# We do not use C#. We do not mess with partial clone either.
RUN git config --global url."https://github.com/".insteadOf "git@github.com:"

WORKDIR /build

COPY Cargo.toml Cargo.lock* build.rs ./
RUN mkdir -p src && echo 'fn main() {}' > src/main.rs && cargo build --release 2>/dev/null || true
RUN rm -rf src
COPY . .
RUN cargo build --release

# ---------------------------------------------------------------------------
# Stage 4: runtime image
# ---------------------------------------------------------------------------
#
# Distroless rather than debian-slim: nothing here links against anything beyond glibc and
# libgcc, and every package a full base drags in (perl, curl/libssh2, wget, gnupg, apt) is
# CVE surface for code this service never runs.
FROM gcr.io/distroless/cc-debian13

COPY --from=builder /build/target/release/grok_connect_adbc /usr/local/bin/

# ClickHouse ADBC driver, built from source in the clickhouse-driver stage
COPY --from=clickhouse-driver /src/target/release/libadbc_clickhouse.so /usr/local/lib/
ENV ADBC_CLICKHOUSE_DRIVER=/usr/local/lib/libadbc_clickhouse.so

# BigQuery and Snowflake ADBC drivers, built in the adbc-go-drivers stage. BigQuery is
# pointed at by full path; Snowflake is loaded by name, hence LD_LIBRARY_PATH — there is
# no ldconfig cache here.
COPY --from=adbc-go-drivers /drivers/libadbc_driver_bigquery.so /usr/local/lib/
COPY --from=adbc-go-drivers /drivers/libadbc_driver_snowflake.so /usr/local/lib/
ENV ADBC_BIGQUERY_DRIVER=/usr/local/lib/libadbc_driver_bigquery.so
ENV LD_LIBRARY_PATH=/usr/local/lib

EXPOSE 1234

HEALTHCHECK --interval=10s --timeout=3s --start-period=5s \
    CMD ["grok_connect_adbc", "--health"]

ENTRYPOINT ["grok_connect_adbc"]
