-- Your SQL goes here
ALTER TABLE `uniprot_entries` DROP COLUMN `mass`;
ALTER TABLE `uniprot_entries` DROP COLUMN `seq_length`;
ALTER TABLE `uniprot_entries` DROP COLUMN `family`;
ALTER TABLE `uniprot_entries` ADD COLUMN `mass` INTEGER NOT NULL;
ALTER TABLE `uniprot_entries` ADD COLUMN `seq_length` INTEGER NOT NULL;
ALTER TABLE `uniprot_entries` ADD COLUMN `family` TEXT NOT NULL;


