fn dydx(k: f64, x: f64, y: f64) -> f64 {
    // -4.0 * k * x.powf(0.0) * y.powf(3.0)
    k * x.powf(0.0) * y
}

fn rk(k: f64, mut x0: f64, x: f64, y0: f64, n: f64) -> f64 { // x0->t0 ; x->t ; y0->Ca0
    let h: f64 = (x-x0) / n; // step size
    let mut y = y0;

    for i in 0..n as u64 {
        let k1 = h * dydx(k, x0, y);
        let k2 = h * dydx(k, x0 + h/2.0, y + k1/2.0);
        let k3 = h * dydx(k, x0 + h/2.0, y + k2/2.0);
        let k4 = h * dydx(k, x0 + h, y + k3); 

        x0 = x0 + h;
        y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    }
    y
}

fn perc_err(x: f64, x_aprox: f64) -> f64 {
    (1.0 - x_aprox/x) * 100.0 // (x - xa) / x
}

fn main() {
    let a = rk(1.0, 0.0, 1.0, 1.0, 1000.0);
    let b = 1.0_f64.exp(); // exp 6
    println!("{}", a);
    println!("{}", b);
    println!("{}", perc_err(b, a));
}