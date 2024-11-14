-- This file should undo anything in `up.sql`
ALTER TABLE `uniprot_entries` DROP COLUMN `mass`;
ALTER TABLE `uniprot_entries` DROP COLUMN `seq_length`;
ALTER TABLE `uniprot_entries` DROP COLUMN `family`;
ALTER TABLE `uniprot_entries` ADD COLUMN `mass` INTEGER;
ALTER TABLE `uniprot_entries` ADD COLUMN `seq_length` INTEGER;
ALTER TABLE `uniprot_entries` ADD COLUMN `family` TEXT;


