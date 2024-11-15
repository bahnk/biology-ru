// @generated automatically by Diesel CLI.

diesel::table! {
    belongs_to_uniprot_sequence_similarity_family (entry, family) {
        entry -> Text,
        family -> Text,
    }
}

diesel::table! {
    uniprot_entries (accession_number) {
        accession_number -> Text,
        entry_name -> Text,
        mass -> Nullable<Integer>,
        seq_length -> Nullable<Integer>,
    }
}

diesel::table! {
    uniprot_sequence_similarity_families (name) {
        name -> Text,
    }
}

diesel::joinable!(belongs_to_uniprot_sequence_similarity_family -> uniprot_entries (entry));
diesel::joinable!(belongs_to_uniprot_sequence_similarity_family -> uniprot_sequence_similarity_families (family));

diesel::allow_tables_to_appear_in_same_query!(
    belongs_to_uniprot_sequence_similarity_family,
    uniprot_entries,
    uniprot_sequence_similarity_families,
);
