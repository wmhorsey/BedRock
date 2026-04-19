use std::collections::HashMap;

use nalgebra::Vector3;

use crate::physics3d::Particle3d;

#[derive(Default)]
pub struct SpatialHashGrid3d {
    cell_size: f32,
    cells: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl SpatialHashGrid3d {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
        }
    }

    pub fn rebuild(&mut self, particles: &[Particle3d]) {
        self.cells.clear();
        for (index, particle) in particles.iter().enumerate() {
            let key = self.cell_key(particle.position);
            self.cells.entry(key).or_default().push(index);
        }
    }

    pub fn query_neighbors(&self, position: Vector3<f32>, out: &mut Vec<usize>) {
        out.clear();
        let (base_x, base_y, base_z) = self.cell_key(position);

        for offset_z in -1..=1 {
            for offset_y in -1..=1 {
                for offset_x in -1..=1 {
                    let key = (base_x + offset_x, base_y + offset_y, base_z + offset_z);
                    if let Some(indices) = self.cells.get(&key) {
                        out.extend(indices.iter().copied());
                    }
                }
            }
        }
    }

    fn cell_key(&self, position: Vector3<f32>) -> (i32, i32, i32) {
        (
            (position.x / self.cell_size).floor() as i32,
            (position.y / self.cell_size).floor() as i32,
            (position.z / self.cell_size).floor() as i32,
        )
    }
}