use diesel::prelude::*;

//#[derive(Queryable, Selectable)]
//#[diesel(table_name = crate::schema::uniprot_entries)]
//#[diesel(check_for_backend(diesel::sqlite::Sqlite))]
pub struct SimilarEntry {
    pub family: String,
    pub entry_name: String,
    pub accession_number: String,
}

//#[derive(Queryable, Selectable)]
//#[diesel(table_name = crate::schema::uniprot_families)]
//#[diesel(check_for_backend(diesel::sqlite::Sqlite))]
pub struct SimilarFamilyEntry {
    pub name: String,
}
