use biology_ru::uniprot::models::*;
use biology_ru::uniprot::similar::get_similar_entries;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let url = "http://www.uniprot.org/docs/similar.txt";

    let entries = get_similar_entries(url)?;

    for entry in entries {
        let (a, b, c) = entry;
        println!("{}", b);
    }

    Ok(())
}
