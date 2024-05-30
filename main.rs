mod coordinate;
mod environment;

use crate::coordinate::{dcm_ned2body_quat, euler2quat};
use crate::environment::{get_std_density, gravity, wind_ned, magnetic_declination};

fn calculate_dynamics(quat: [f64; 4], altitude: f64) -> f64 {
    let dcm = dcm_ned2body_quat(quat);
    let density = get_std_density(altitude);
    let g = gravity(altitude);
    let result = dcm[0][0] * density * g;
    result
}

fn main() {
    let quat = [1.0, 0.0, 0.0, 0.0];
    let altitude = 1000.0;
    let result = calculate_dynamics(quat, altitude);
    println!("計算結果: {}", result);
}
