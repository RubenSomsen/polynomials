// Learning project to better understand polynomials
// and how they relate to threshold signature schemes

// Lagrange interpolation is fully functional

// TODO: needs refactoring into modular arithmetic
// TODO: polynomial interpolation (get poly from points)
// TODO: study up on lagrange coefficients

fn main() {
    println!("NOW RUNNING");
    
    use std::time::Instant;
    let now = Instant::now();

    calc_polynomial();

    let elapsed = now.elapsed();
    println!("PERFORMANCE: {:?}", elapsed);
}

fn calc_polynomial() {
    let p1 = Polynomial::new(vec![7.,3.]); // degree 1
    let p2 = Polynomial::new_random(2); // degree 2
    let p3 = Polynomial::new(vec![9.,3.,2.,5.]); // degree 3
    
    // Get y values from polynomial
    let y1 = p1.get_y(1.);
    let y2 = p1.get_y(2.);
    let input = vec![y1, y2];
    assert_eq!(input, p1.get_y_values());

    // Given degree+1 y coordinates, calculate any other point on the curve
    let target_x = 4.;
    let secret = Polynomial::lagrange_interpolate(&input, target_x);
    println!("Secret: {} {}", secret, p1.get_y(4.));
    assert_eq!(secret, p1.get_y(target_x));

    p1.print();
    let v = Polynomial::lagrange_interpolate(&p1.get_y_values(), 2.);
    println!("{} {} {}", p3.get_y(6.), p1.get_y(6.), Polynomial::gcd(4,18));
    return
}

struct Polynomial {
    coefficients: Vec<f64>
}

impl Polynomial {
    fn new(coefficients: Vec<f64>) -> Polynomial {
        return Polynomial { coefficients }
    }

    fn new_random(degree: usize) -> Polynomial {
        use rand::{distributions::Uniform, Rng};
        let mut rng = rand::thread_rng();
        let range = Uniform::new(0., f32::MAX as f64);
        return Polynomial::new((0..degree+1).map(|_| rng.sample(&range)).collect())
    }

    fn get_y(&self, x: f64) -> f64 {
        let mut y = 0.;
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            y += coefficient*x.powf(i as f64);
        }
        return y
    }

    fn get_y_values(&self) -> Vec<f64> {
        let mut y_values = vec![];
        for i in 0..self.coefficients.len() {
            y_values.push(self.get_y(i as f64+1.));
        }
        return y_values
    }

    // Given degree+1 y coordinates, calculate any other y coordinate on a given x coordinate
    fn lagrange_interpolate(y_values: &Vec<f64>, target_x: f64) -> f64 {
        let mut y = 0.;
        for x in 1..1+y_values.len() { // assumes x == index+1 (TODO: clean up)
            let (mut dividend, mut divisor) = (1., 1.);
            for xx in 1..1+y_values.len() {
                if x==xx { continue }
                dividend *= target_x-xx as f64;
                divisor *= x as f64-xx as f64;
            }
            y += (y_values[x as usize-1]*dividend)/divisor;
        }
        return y
    }

    // Given degree+1 y coordinates, recreate the polynomial
    // TODO: finish function
    fn interpolate_polynomial(y_values: &Vec<f64>) -> Polynomial {
        let mut y_values = y_values.clone();
        let c0 = Polynomial::lagrange_interpolate(&y_values, 0.);
        let mut coefficients = vec![];
        y_values[0] -= c0;
        //let subtractors = vec![coefficients[0]-c0];
        for x in 2..1+y_values.len() { // assumes x == index+1
            let y = y_values[x-1];

            // TODO
            // It's actually quite a few steps, basically gaussian elimination
            // maybe OK, iter down and then back up
        }
        coefficients.push(c0);
        coefficients.reverse();
        return Polynomial::new(coefficients);
    }

    // Greatest common divisor  
    fn gcd(a: u64, b: u64) -> u64 { // a is smaller
        let (a, b) = if a > b { (b, a) } else { (a, b) };
        if a == 0 { return b }
        let r = b % a; // remainder
        return Polynomial::gcd(r, a)
    }

    // Prints the polynomial
    fn print(&self) {
        for (i, coefficient) in self.coefficients.iter().enumerate() {
            if i == 0 { print!("{}", coefficient) }
            else { print!(" + {}*x^{}", coefficient, i) }
        }
        println!();
    }
}