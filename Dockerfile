FROM rust:1.94-trixie AS build

# Debian trixie ships htslib 1.21; the core needs 1.18 or newer.
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    libhts-dev \
    pkg-config \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cargo build --release --locked -p mito-cli -p mito-server

FROM debian:trixie-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libhts3t64 \
  && rm -rf /var/lib/apt/lists/* \
  && useradd --system --uid 10001 --home-dir /app --shell /usr/sbin/nologin mito \
  && mkdir -p /app/tmp \
  && chown mito:mito /app/tmp

WORKDIR /app
COPY --from=build /src/target/release/mito-cli /usr/local/bin/mito-cli
COPY --from=build /src/target/release/mito-server /usr/local/bin/mito-server
COPY --from=build /src/core/data /app/core/data

# The binaries embed build-tree resource paths; point them at this copy.
ENV MITO_DATA_DIR=/app/core/data
ENV MITO_SERVER_ADDR=0.0.0.0:8080
USER mito
EXPOSE 8080
CMD ["mito-server"]
