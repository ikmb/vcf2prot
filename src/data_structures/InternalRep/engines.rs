use std::str::FromStr;

/// An enumerator describing possible cases to generate an execution engine 
/// Currently, three engines are supported, a single-threaded engine, a multi-threaded engine and a GPU execution engine 
/// the enumerator implements the FromStr trait and derives the Debug and Clone trait 
/// ```rust
/// let engine=Engine::from_str("st"); 
/// match engine
/// {
///     Engine::ST=>println!("An Engine instance of type: {:#?}",engine),
///     _=>println!("This line will not be printed")
/// }
/// ´´´
#[derive(Debug,Clone)]
pub enum Engine{ST,MT,GPU}

impl FromStr for Engine 
{
    type Err=String;
    fn from_str(eninge_name:&str)->Result<Engine,String>
    {
        match eninge_name
        {
            "st"  | "ST" =>Ok(Engine::ST), 
            "mt"  | "MT" =>Ok(Engine::MT),
            "gpu" | "GPU"=>Ok(Engine::GPU),
            _=>Err(format!("{} is not a supported engine",eninge_name))
        }
    }
}