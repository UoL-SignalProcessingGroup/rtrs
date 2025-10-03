
// math helpers
// #[inline] pub fn sign_change(a: f32, b: f32) -> bool { (a<=0.0 && b>0.0) || (a>=0.0 && b<0.0) }
#[inline] pub fn dot(a: &[f32;3], b: &[f32;3]) -> f32 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
#[inline] pub fn norm(a: &[f32;3]) -> f32 { dot(a,a).sqrt() }
#[inline] pub fn sub(a: &[f32;3], b: &[f32;3]) -> [f32;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
// #[inline] pub fn sub2(a: &[f32;2], b: &[f32;2]) -> [f32;2] { [a[0]-b[0], a[1]-b[1]] }
// #[inline] pub fn muls(a: &[f32;3], s: f32) -> [f32;3] { [a[0]*s, a[1]*s, a[2]*s] }
#[inline] pub fn cross_product(a: &[f32;3], b: &[f32;3]) -> [f32;3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
// #[inline] pub fn add(a: &[f32;3], b: &[f32;3]) -> [f32;3] { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
// #[inline] pub fn cross(a: &[f32;3], b: &[f32;3]) -> [f32;3] {
//     [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
// }