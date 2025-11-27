FROM node:18-alpine AS mcp-build

RUN apk add --no-cache git

WORKDIR /app

RUN git clone https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git
WORKDIR /app/ChEMBL-MCP-Server

RUN npm ci
RUN npm run build

FROM node:18-alpine

WORKDIR /app

RUN apk add --no-cache curl unzip
RUN curl -sLo ngrok.zip https://bin.equinox.io/c/bNyj1mQVY4c/ngrok-v3-stable-linux-amd64.tgz \
    && tar -xzf ngrok.zip -C /usr/local/bin \
    && rm ngrok.zip

COPY --from=mcp-build /app/ChEMBL-MCP-Server ./chembl-mcp-server

COPY package*.json ./
COPY tsconfig.json ./
RUN npm i

COPY src/ ./src/
RUN npm run build

RUN addgroup -g 1001 -S appuser && adduser -S appuser -u 1001
RUN chown -R appuser:appuser /app
USER appuser

EXPOSE 3333

ENV NGROK_AUTH_TOKEN="360eXoKoi6t9vosMqoWZPGMicQw_3kP2wVax8KdYfdJXBr6S6"

ENTRYPOINT sh -c "ngrok authtoken $NGROK_AUTH_TOKEN && ngrok http 3333 --log=stdout & node build/index.js"
