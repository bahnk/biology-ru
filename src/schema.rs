// @generated automatically by Diesel CLI.

diesel::table! {
    uniprot_entries (accession_number) {
        accession_number -> Text,
        entry_name -> Text,
        mass -> Integer,
        seq_length -> Integer,
        family -> Text,
    }
}

diesel::table! {
    uniprot_families (name) {
        name -> Text,
    }
}

diesel::joinable!(uniprot_entries -> uniprot_families (family));

diesel::allow_tables_to_appear_in_same_query!(
    uniprot_entries,
    uniprot_families,
);
