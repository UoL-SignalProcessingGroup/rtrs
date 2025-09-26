
// math helpers
// #[inline] pub fn sign_change(a: f64, b: f64) -> bool { (a<=0.0 && b>0.0) || (a>=0.0 && b<0.0) }
#[inline] pub fn dot(a: &[f64;3], b: &[f64;3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
#[inline] pub fn norm(a: &[f64;3]) -> f64 { dot(a,a).sqrt() }
#[inline] pub fn sub(a: &[f64;3], b: &[f64;3]) -> [f64;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
// #[inline] pub fn sub2(a: &[f64;2], b: &[f64;2]) -> [f64;2] { [a[0]-b[0], a[1]-b[1]] }
// #[inline] pub fn muls(a: &[f64;3], s: f64) -> [f64;3] { [a[0]*s, a[1]*s, a[2]*s] }
#[inline] pub fn cross_product(a: &[f64;3], b: &[f64;3]) -> [f64;3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
// #[inline] pub fn add(a: &[f64;3], b: &[f64;3]) -> [f64;3] { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
// #[inline] pub fn cross(a: &[f64;3], b: &[f64;3]) -> [f64;3] {
//     [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
// }