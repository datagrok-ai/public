-- Table for storing canonical SMILES with auto-generated VIDs
CREATE TABLE hitdesign.vid_dictionary (
    id SERIAL PRIMARY KEY,
    vid VARCHAR(15) GENERATED ALWAYS AS ('V' || LPAD(id::TEXT, 6, '0')) STORED UNIQUE,
    mh_string TEXT NOT NULL UNIQUE
);

-- Table for tracking which VIDs are used in which campaigns
CREATE TABLE hitdesign.campaign_vids (
    app_name VARCHAR(100) NOT NULL,
    vid VARCHAR(15) NOT NULL REFERENCES hitdesign.vid_dictionary(vid),
    campaign_id VARCHAR(255) NOT NULL,
    created_by TEXT,
    PRIMARY KEY (app_name, vid, campaign_id)
);

