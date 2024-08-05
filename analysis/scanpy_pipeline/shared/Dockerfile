ARG PIXI_VERSION=0.26.1
ARG BASE_IMAGE=debian:bookworm-slim

FROM --platform=$TARGETPLATFORM ubuntu:24.04 AS builder
# need to specify the ARG again to make it available in this stage
ARG PIXI_VERSION
RUN apt-get update && apt-get install -y curl
# download the musl build since the gnu build is not available on aarch64
RUN curl -Ls \
    "https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl" \
    -o /pixi && chmod +x /pixi
RUN /pixi --version

FROM --platform=$TARGETPLATFORM $BASE_IMAGE
COPY --from=builder --chown=root:root --chmod=0555 /pixi /usr/local/bin/pixi
# copy source code, pixi.toml and pixi.lock to the container

# Set working directory
WORKDIR /app
RUN mkdir /app/host

# Set permissions if needed
RUN chmod 770 /app/host

# Copy pixi.toml file
COPY pixi.toml .

# Install dependencies using pixi
RUN pixi install

# Copy your Python script
COPY sc_conversion.R .
COPY setupenv.R .
COPY assets/GenomeInfoDbData_1.2.11.tar.gz .

# Copy the entrypoint script
COPY entrypoint.sh .

RUN pixi run "R CMD INSTALL GenomeInfoDbData_1.2.11.tar.gz"
RUN pixi run "Rscript setupenv.R"

# Make the entrypoint script executable
RUN chmod +x entrypoint.sh

# Set the entrypoint
ENTRYPOINT ["./entrypoint.sh"]

