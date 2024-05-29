extern crate nalgebra as na;

use std::f64::consts::PI;

fn dcm_wind2body(alpha: f64, beta: f64) -> na::Matrix3<f64> {
    let dcm_0 = na::Vector3::new(alpha.cos() * beta.cos(), alpha.cos() * beta.sin(), -alpha.sin());
    let dcm_1 = na::Vector3::new(-beta.sin(), beta.cos(), 0.0);
    let dcm_2 = na::Vector3::new(alpha.sin() * beta.cos(), alpha.sin() * beta.sin(), alpha.cos());
    na::Matrix3::from_columns(&[dcm_0, dcm_1, dcm_2])
}

fn dcm_ned2body_euler(azimuth: f64, elevation: f64, roll: f64) -> na::Matrix3<f64> {
    let azi = azimuth;
    let elv = elevation;
    let rol = roll;
    let dcm_0 = na::Vector3::new(azi.cos() * elv.cos(), azi.sin() * elv.cos(), -elv.sin());
    let dcm_1 = na::Vector3::new(
        -azi.sin() * rol.cos() + azi.cos() * elv.sin() * rol.sin(),
        azi.cos() * rol.cos() + azi.sin() * elv.sin() * rol.sin(),
        elv.cos() * rol.sin(),
    );
    let dcm_2 = na::Vector3::new(
        azi.sin() * rol.sin() + azi.cos() * elv.sin() * rol.cos(),
        -azi.cos() * rol.sin() + azi.sin() * elv.sin() * rol.cos(),
        elv.cos() * rol.cos(),
    );
    na::Matrix3::from_columns(&[dcm_0, dcm_1, dcm_2])
}

fn quat_normalize(quat: na::Quaternion<f64>) -> na::Quaternion<f64> {
    quat.normalize()
}

fn dcm_ned2body_quat(quat: na::Quaternion<f64>) -> na::Matrix3<f64> {
    let q0 = quat.w;
    let q1 = quat.i;
    let q2 = quat.j;
    let q3 = quat.k;

    let dcm_0 = na::Vector3::new(q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3, 2.0 * (q0 * q1 + q2 * q3), 2.0 * (q0 * q2 - q1 * q3));
    let dcm_1 = na::Vector3::new(2.0 * (q0 * q1 - q2 * q3), q1 * q1 - q0 * q0 - q2 * q2 + q3 * q3, 2.0 * (q1 * q2 + q0 * q3));
    let dcm_2 = na::Vector3::new(2.0 * (q0 * q2 + q1 * q3), 2.0 * (q1 * q2 - q0 * q3), q2 * q2 - q0 * q0 - q1 * q1 + q3 * q3);
    na::Matrix3::from_columns(&[dcm_0, dcm_1, dcm_2])
}

fn euler2quat(azimuth: f64, elevation: f64, roll: f64) -> na::Quaternion<f64> {
    let azimuth = azimuth.to_radians();
    let elevation = elevation.to_radians();
    let roll = roll.to_radians();

    let dcm = dcm_ned2body_euler(azimuth, elevation, roll);
    let q0 = 0.5 * (1.0 + dcm[(0, 0)] - dcm[(1, 1)] - dcm[(2, 2)]).sqrt();
    let q1 = 0.5 * (1.0 - dcm[(0, 0)] + dcm[(1, 1)] - dcm[(2, 2)]).sqrt();
    let q2 = 0.5 * (1.0 - dcm[(0, 0)] - dcm[(1, 1)] + dcm[(2, 2)]).sqrt();
    let q3 = 0.5 * (1.0 + dcm[(0, 0)] + dcm[(1, 1)] + dcm[(2, 2)]).sqrt();

    let quat = na::Quaternion::new(q0, q1, q2, q3);
    quat.normalize()
}

fn quat2euler(dcm_ned2body: na::Matrix3<f64>) -> (f64, f64, f64) {
    let dcm = dcm_ned2body;
    let azimuth = (dcm[(0, 1)]).atan2(dcm[(0, 0)]).to_degrees();
    let elevation = (-dcm[(0, 2)]).asin().to_degrees();
    let roll = (dcm[(1, 2)]).atan2(dcm[(2, 2)]).to_degrees();

    (azimuth, elevation, roll)
}

fn dcm_eci2ecef(t_sec: f64) -> na::Matrix3<f64> {
    let omega_e = 7.2921151e-5; // [rad/s] for WGS84
    let xi = omega_e * t_sec;
    let dcm_0 = na::Vector3::new(xi.cos(), xi.sin(), 0.0);
    let dcm_1 = na::Vector3::new(-xi.sin(), xi.cos(), 0.0);
    let dcm_2 = na::Vector3::new(0.0, 0.0, 1.0);
    na::Matrix3::from_columns(&[dcm_0, dcm_1, dcm_2])
}

fn dcm_ecef2ned(pos0_llh: na::Vector3<f64>) -> na::Matrix3<f64> {
    let lat = pos0_llh[0].to_radians();
    let lon = pos0_llh[1].to_radians();
    let dcm_0 = na::Vector3::new(-lat.sin() * lon.cos(), -lat.sin() * lon.sin(), lat.cos());
    let dcm_1 = na::Vector3::new(-lon.sin(), lon.cos(), 0.0);
    let dcm_2 = na::Vector3::new(-lat.cos() * lon.cos(), -lat.cos() * lon.sin(), -lat.sin());
    na::Matrix3::from_columns(&[dcm_0, dcm_1, dcm_2])
}

fn vel_eci2ecef(vel_eci: na::Vector3<f64>, dcm_eci2ecef: na::Matrix3<f64>, pos_eci: na::Vector3<f64>) -> na::Vector3<f64> {
    let omega_e = 7.2921151e-5; // [rad/s] for WGS84
    let omega_eci2ecef = na::Matrix3::new(0.0, -omega_e, 0.0, omega_e, 0.0, 0.0, 0.0, 0.0, 0.0);
    dcm_eci2ecef * vel_eci - omega_eci2ecef * pos_eci
}

fn vel_ecef2eci(vel_ecef: na::Vector3<f64>, dcm_eci2ecef: na::Matrix3<f64>, pos_eci: na::Vector3<f64>) -> na::Vector3<f64> {
    let omega_e = 7.2921151e-5; // [rad/s] for WGS84
    let omega_eci2ecef = na::Matrix3::new(0.0, -omega_e, 0.0, omega_e, 0.0, 0.0, 0.0, 0.0, 0.0);
    dcm_eci2ecef.transpose() * vel_ecef + omega_eci2ecef * pos_eci
}
