use std::str::FromStr;

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