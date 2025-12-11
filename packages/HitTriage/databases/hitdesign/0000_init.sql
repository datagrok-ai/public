-- Simple lock table
CREATE TABLE hitdesign.campaign_locks (
    app_name VARCHAR(100) NOT NULL,
    campaign_id VARCHAR(255) NOT NULL,
    expires_at TIMESTAMP WITH TIME ZONE NOT NULL,
    --optional: user who holds the lock
    locked_by TEXT,
    PRIMARY KEY (app_name, campaign_id)
);

CREATE TABLE hitdesign.update_logs (
    app_name VARCHAR(100) NOT NULL,
    campaign_id VARCHAR(255) NOT NULL,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (app_name, campaign_id)
);

GRANT ALL PRIVILEGES ON SCHEMA hitdesign TO CURRENT_USER;
GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA hitdesign TO :LOGIN;
GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA hitdesign TO :LOGIN;