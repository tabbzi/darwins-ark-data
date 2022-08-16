## Repository

### darwins-ark-data

Extracting, cleaning, and preparing raw data tables from Darwin's Ark database for downstream analysis

#### Secrets

The following secret variables are used in this repository:

- `DB_CERTIFICATE`: cerificate for database access
- `DB_USERNAME`: username for database access
- `DB_PASSWORD`: password for database access

## Project

### Darwin's Ark

[Darwin's Ark](darwinsark.org) is a community science initiative and platform for collecting genetic, environmental, health, and behavioral data from companion and working animals. Darwin's Dogs, its largest study, focuses on the genomics of health and behavior in companion dogs.

### Database

The project maintains a SQL database for storing data. Access to the database is private and requires a certificate, username, and password via the `mysql` command-line tool or a staff account to access via phpMyAdmin. The command-line SQL commands used to extract tables in `/raw/` are provided in `/sql/`. The directory `/raw/YYYYMMDD/` contains examples of raw data tables exported from the database using these commands are given for staff owners and dogs, including test dogs (`fake==1`).