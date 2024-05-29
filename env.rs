extern crate nalgebra as na;

use std::f64::consts::PI;

fn std_atmo(altitude: f64) -> (f64, f64, f64, f64) {
    let r = 287.1;
    let gamma = 1.4;
    let re = 6378.137e3;
    let g0 = 9.80665;

    let h_list = [0.0, 11.0e3, 20.0e3, 32.0e3, 47.0e3, 51.0e3, 71.0e3, 84.852e3];
    let tg_list = [-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3, 0.0];
    let t_list = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.946];
    let p_list = [101325.0, 22632.0, 5474.9, 868.02, 110.91, 66.939, 3.9564, 0.3734];

    let h = altitude * re / (re + altitude);

    let mut k = 0;
    for i in 0..8 {
        if h < h_list[i] {
            k = if i == 0 { 0 } else { i - 1 };
            break;
        } else if h >= h_list[7] {
            k = 7;
            break;
        }
    }

    let temperature = t_list[k] + tg_list[k] * (h - h_list[k]);
    let pressure = if tg_list[k] == 0.0 {
        p_list[k] * f64::exp(g0 / r * (h_list[k] - h) / t_list[k])
    } else {
        p_list[k] * f64::powf(t_list[k] / temperature, g0 / r / tg_list[k])
    };
    let density = pressure / (r * temperature);
    let sound_speed = f64::sqrt(gamma * r * temperature);

    (temperature, pressure, density, sound_speed)
}

fn get_std_temp(altitude: f64) -> f64 {
    std_atmo(altitude).0
}

fn get_std_temp_array(altitude_array: &[f64]) -> Vec<f64> {
    altitude_array.iter().map(|&alt| get_std_temp(alt)).collect()
}

fn get_std_press(altitude: f64) -> f64 {
    std_atmo(altitude).1
}

fn get_std_press_array(altitude_array: &[f64]) -> Vec<f64> {
    altitude_array.iter().map(|&alt| get_std_press(alt)).collect()
}

fn get_std_density(altitude: f64) -> f64 {
    std_atmo(altitude).2
}

fn get_std_density_array(altitude_array: &[f64]) -> Vec<f64> {
    altitude_array.iter().map(|&alt| get_std_density(alt)).collect()
}

fn get_std_soundspeed(altitude: f64) -> f64 {
    std_atmo(altitude).3
}

fn get_std_soundspeed_array(altitude_array: &[f64]) -> Vec<f64> {
    altitude_array.iter().map(|&alt| get_std_soundspeed(alt)).collect()
}

fn gravity(altitude: f64) -> f64 {
    let re = 6378.137e3;
    let g0 = 9.80665;
    let alt = if altitude < 0.0 { 0.0 } else { altitude };
    g0 * f64::powf(re / (re + alt), 2.0)
}

fn wind_ned(wind_speed: f64, wind_direction: f64) -> na::Vector3<f64> {
    let wind_speed = -wind_speed;
    let wind_direction_rad = wind_direction.to_radians();
    na::Vector3::new(
        wind_speed * wind_direction_rad.cos(),
        wind_speed * wind_direction_rad.sin(),
        0.0,
    )
}

fn magnetic_declination(lat: f64, lon: f64) -> f64 {
    let delta_lat = lat - 37.0;
    let delta_lon = lon - 138.0;
    (7.0 + 57.201 / 60.0)
        + (18.750 / 60.0) * delta_lat
        - (6.761 / 60.0) * delta_lon
        - (0.059 / 60.0) * delta_lat.powi(2)
        - (0.014 / 60.0) * delta_lat * delta_lon
        - (0.579 / 60.0) * delta_lon.powi(2)
}

fn main() {
    let altitude = 1000.0;
    let (temp, press, dens, sound_speed) = std_atmo(altitude);
    println!("Temperature: {} K", temp);
    println!("Pressure: {} Pa", press);
    println!("Density: {} kg/m^3", dens);
    println!("Sound Speed: {} m/s", sound_speed);

    let altitudes = vec![0.0, 1000.0, 2000.0];
    let temps = get_std_temp_array(&altitudes);
    println!("Temperatures: {:?}", temps);

    let winds = wind_ned(10.0, 45.0);
    println!("Wind NED: {:?}", winds);

    let mag_dec = magnetic_declination(35.0, 139.0);
    println!("Magnetic Declination: {} degrees", mag_dec);
}
