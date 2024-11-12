CREATE TABLE uniprot_entries (
  -- Entry, accession
  accession_number VARCHAR(20) NOT NULL PRIMARY KEY,

  -- Entry Name, id
  entry_name VARCHAR(100) NOT NULL,

  -- Mass, mass
  mass INTEGER NOT NULL,

  -- Length, length
  seq_length INTEGER NOT NULL,

  -- Similarity family
  family VARCHAR NOT NULL,
  FOREIGN KEY (family) REFERENCES uniprot_families(name)
)
