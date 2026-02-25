
// math helpers
#[inline] pub fn dot(a: &[f32;3], b: &[f32;3]) -> f32 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
#[inline] pub fn sub(a: &[f32;3], b: &[f32;3]) -> [f32;3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }