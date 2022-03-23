@echo Assure running the Docker daemon (Docker Desktop)
docker compose --project-name datagrok -f="localhost.docker-compose.yaml" --profile all up