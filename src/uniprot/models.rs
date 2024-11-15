use diesel::prelude::*;

#[derive(Queryable, Selectable, Insertable, Clone)]
#[diesel(table_name = crate::schema::uniprot_entries)]
#[diesel(check_for_backend(diesel::sqlite::Sqlite))]
pub struct UniprotEntry {
    pub entry_name: String,
    pub accession_number: String,
    pub mass: Option<i32>,
    pub seq_length: Option<i32>,
}

#[derive(Queryable, Selectable, Insertable, Clone)]
#[diesel(table_name = crate::schema::uniprot_sequence_similarity_families)]
#[diesel(check_for_backend(diesel::sqlite::Sqlite))]
pub struct UniprotFamily {
    pub name: String,
}

#[derive(Queryable, Selectable, Insertable)]
#[diesel(table_name = crate::schema::belongs_to_uniprot_sequence_similarity_family)]
#[diesel(check_for_backend(diesel::sqlite::Sqlite))]
pub struct BelongsToFamily {
    pub entry: String,
    pub family: String,
}
