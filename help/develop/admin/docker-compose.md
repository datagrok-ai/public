
<!-- TITLE: Deploy Datagrok using Docker Compose -->
<!-- SUBTITLE: -->

# Deploy Datagrok using Docker Compose

This document contains instructions to deploy Datagrok on a regular machine without cloud-based hosting.

Create `docker-compose.yaml` in an empty folder
```yaml
version: "3"
services:
  db:
    image: postgres
    environment:
      POSTGRES_USER: postgres
      POSTGRES_PASSWORD: postgres
    networks:
      datagrok:
        aliases:
          - database
    volumes:
      - datagrok_db:/var/lib/postgresql/data
  datagrok:
    image: datagrok/datagrok:1.0.82-75b821b
    environment:
      GROK_PARAMETERS: "{\"deployDemo\": false, \"dbServer\": \"database\", \"db\": \"datagrok\", \"dbAdminLogin\": \"postgres\", \"dbAdminPassword\": \"postgres\", \"dbLogin\": \"dg\", \"dbPassword\": \"dg\"}"
    ports:
      - "8080:8080"
    networks:
      datagrok:
        aliases:
          - datagrok
    volumes:
      - datagrok_data:/home/grok/data
      - datagrok_cfg:/home/grok/cfg
  cvm:
    image: datagrok/cvm:1.0.82-75b821b
    environment:
      GROK_COMPUTE_NUM_CORES: 4
    ports:
      - "5005:5005"
      - "8090:8090"
      - "54321:54321"
    networks:
      datagrok:
        aliases:
          - cvm
volumes: 
  datagrok_db:
  datagrok_data:
  datagrok_cfg:
networks:
  datagrok:
```

Run `docker-compose up` command

Datagrok will start and automatically deploy a new database. After it deployed first time you can shut it down using `Ctrl+C`.
All data will be saved in persistent storage.

See also:

* [Docker Compose](https://docs.docker.com/compose/)