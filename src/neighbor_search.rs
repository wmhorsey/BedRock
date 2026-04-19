use std::collections::HashMap;

use nalgebra::Vector2;

use crate::physics::Particle;

#[derive(Default)]
pub struct SpatialHashGrid {
    cell_size: f32,
    cells: HashMap<(i32, i32), Vec<usize>>,
}

impl SpatialHashGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cells: HashMap::new(),
        }
    }

    pub fn rebuild(&mut self, particles: &[Particle]) {
        self.cells.clear();
        for (index, particle) in particles.iter().enumerate() {
            let key = self.cell_key(particle.position);
            self.cells.entry(key).or_default().push(index);
        }
    }

    pub fn query_neighbors(&self, position: Vector2<f32>, out: &mut Vec<usize>) {
        out.clear();
        let (base_x, base_y) = self.cell_key(position);

        for offset_y in -1..=1 {
            for offset_x in -1..=1 {
                let key = (base_x + offset_x, base_y + offset_y);
                if let Some(indices) = self.cells.get(&key) {
                    out.extend(indices.iter().copied());
                }
            }
        }
    }

    fn cell_key(&self, position: Vector2<f32>) -> (i32, i32) {
        (
            (position.x / self.cell_size).floor() as i32,
            (position.y / self.cell_size).floor() as i32,
        )
    }
}