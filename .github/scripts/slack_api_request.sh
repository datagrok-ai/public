#!/bin/bash

#set -x

# Function to make a Slack API request and handle rate limiting (status 429)
slack_api_request() {
  local url=$1   # Slack API URL
  local method=$2   # HTTP method (GET, POST, etc.)
  local token=$3   # Slack token
  local data=$4   # Optional JSON data for the request

  MAX_RETRIES=5
  RETRY_COUNT=0
  SUCCESS=false

  while [ "$SUCCESS" = false ] && [ "$RETRY_COUNT" -lt "$MAX_RETRIES" ]; do
#    echo "Making request to Slack API: $url (Attempt #$((RETRY_COUNT+1)))"

    # Send API request
    if [ "$method" = "POST" ]; then
      RESPONSE=$(curl -s -w "%{http_code}" -X POST "$url" \
        -H "Authorization: Bearer $token" \
        -H "Content-Type: application/json; charset=utf-8" \
        --data "$data")
    else
      RESPONSE=$(curl -s -w "%{http_code}" -X "$method" "$url$data" \
        -H "Authorization: Bearer $token" \
        -H "Content-Type: application/json; charset=utf-8")
    fi

    # Extract HTTP status code
    HTTP_STATUS="${RESPONSE: -3}"
    BODY="${RESPONSE%$HTTP_STATUS}"

    if [ "$HTTP_STATUS" -eq 429 ]; then
      # Rate limit hit, extract Retry-After header and wait
#      echo "Rate limit hit..."
      sleep 5
      RETRY_COUNT=$((RETRY_COUNT+1))
    else
      echo "$BODY"
      SUCCESS=true
    fi
  done

  if [ "$SUCCESS" = false ]; then
    echo "Failed to complete the request after $MAX_RETRIES retries due to rate limits."
    exit 1
  fi
}

# Check if parameters are provided correctly
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: slack_api_request.sh <API_URL> <HTTP_METHOD> <SLACK_TOKEN> [<JSON_DATA>]"
  exit 1
fi

# Extract input parameters
SLACK_API_URL="$1"
HTTP_METHOD="$2"
SLACK_TOKEN="$3"
JSON_DATA="$4"

# Call the Slack API request function
slack_api_request "$SLACK_API_URL" "$HTTP_METHOD" "$SLACK_TOKEN" "$JSON_DATA"
