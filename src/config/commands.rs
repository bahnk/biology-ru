use serde::Deserialize;

use super::uaspire;
use super::uniprot;

#[derive(Debug, Deserialize)]
pub enum Commands {
    Uaspire(uaspire::Config),
    Uniprot(uniprot::Config),
}
