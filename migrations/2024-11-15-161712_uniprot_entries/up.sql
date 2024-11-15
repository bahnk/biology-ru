CREATE TABLE uniprot_entries (
  -- Entry, accession
  accession_number VARCHAR(50) NOT NULL PRIMARY KEY,

  -- Entry Name, id
  entry_name VARCHAR(50) NOT NULL,

  -- Mass, mass
  mass INTEGER,

  -- Length, length
  seq_length INTEGER
)
