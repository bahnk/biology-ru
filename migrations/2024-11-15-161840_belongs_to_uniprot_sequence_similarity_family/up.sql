CREATE TABLE belongs_to_uniprot_sequence_similarity_family (
  entry VARCHAR(50) NOT NULL,
  family VARCHAR(300) NOT NULL,
  PRIMARY KEY (entry, family),
  FOREIGN KEY (entry) REFERENCES uniprot_entries(accession_number),
  FOREIGN KEY (family) REFERENCES uniprot_sequence_similarity_families(name)
)
