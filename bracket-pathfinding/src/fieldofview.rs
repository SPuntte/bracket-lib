use std::{collections::HashSet, iter::Map, ops::RangeInclusive};

use num_rational::Rational32;

use bracket_algorithm_traits::prelude::Algorithm2D;
use bracket_geometry::prelude::Point;

enum Cardinal {
    North,
    East,
    South,
    West,
}

struct Quadrant {
    cardinal: Cardinal,
    origin: Point,
}

impl Quadrant {
    pub fn new(cardinal: Cardinal, origin: Point) -> Self {
        Self { cardinal, origin }
    }

    pub fn transform(&self, tile: Tile) -> Point {
        match self.cardinal {
            Cardinal::North => Point::new(self.origin.x + tile.column, self.origin.y - tile.depth),
            Cardinal::East => Point::new(self.origin.x + tile.column, self.origin.y + tile.depth),
            Cardinal::South => Point::new(self.origin.x + tile.depth, self.origin.y + tile.column),
            Cardinal::West => Point::new(self.origin.x - tile.depth, self.origin.y + tile.column),
        }
    }
}

struct Scanline {
    pub depth: i32,
    pub start_slope: Rational32,
    pub end_slope: Rational32,
}

#[derive(Clone, Copy)]
struct Tile {
    pub depth: i32,
    pub column: i32,
}

impl Scanline {
    fn with_integers(depth: i32, start_slope: i32, end_slope: i32) -> Self {
        Self::new(
            depth,
            Rational32::from_integer(start_slope),
            Rational32::from_integer(end_slope),
        )
    }

    fn new(depth: i32, start_slope: Rational32, end_slope: Rational32) -> Self {
        Self {
            depth,
            start_slope,
            end_slope,
        }
    }

    fn tiles(&self) -> Map<RangeInclusive<i32>, impl FnMut(i32) -> Tile> {
        let start_column = round_ties_up(Rational32::from_integer(self.depth) * self.start_slope);
        let end_column = round_ties_down(Rational32::from_integer(self.depth) * self.end_slope);
        let depth = self.depth;
        (start_column..=end_column)
            .into_iter()
            .map(move |column| Tile { depth, column })
    }

    fn next(&self) -> Self {
        Self::new(self.depth + 1, self.start_slope, self.end_slope)
    }
}

fn slope(tile: Tile) -> Rational32 {
    Rational32::new(2 * tile.column - 1, 2 * tile.depth)
}

fn round_ties_up(r: Rational32) -> i32 {
    (r + Rational32::new(1, 2)).floor().to_integer()
}

fn round_ties_down(r: Rational32) -> i32 {
    (r - Rational32::new(1, 2)).ceil().to_integer()
}

fn is_symmetric(scanline: &Scanline, tile: Tile) -> bool {
    let column = Rational32::from_integer(tile.column);
    let depth = Rational32::from_integer(scanline.depth);
    column >= depth * scanline.start_slope && column <= depth * scanline.end_slope
}

struct FovScanner<'a, T: Algorithm2D> {
    radius_2: i32,
    radius_plus_half_2: Rational32,
    quadrant: Quadrant,
    fov_check: &'a T,
    visible_points: &'a mut HashSet<Point>,
}

impl<T: Algorithm2D> FovScanner<'_, T> {
    fn reveal(&mut self, tile: Tile) {
        self.visible_points.insert(self.quadrant.transform(tile));
    }

    fn is_opaque(&mut self, tile: Tile) -> bool {
        let point = self.quadrant.transform(tile);
        if self.fov_check.in_bounds(point) {
            self.fov_check
                .is_opaque(self.fov_check.point2d_to_index(point))
        } else {
            true
        }
    }

    fn scan(&mut self, first_line: Scanline) {
        let mut stack = vec![first_line];
        let mut prev_tile;
        while let Some(mut scanline) = stack.pop() {
            if scanline.depth * scanline.depth > self.radius_2 {
                continue;
            }
            prev_tile = None;
            for tile in scanline.tiles() {
                let tile_point = self.quadrant.transform(tile);
                let dx = tile_point.x - self.quadrant.origin.x;
                let dy = tile_point.y - self.quadrant.origin.y;
                if Rational32::from_integer(dx * dx + dy * dy) <= self.radius_plus_half_2 {
                    if self.is_opaque(tile) || is_symmetric(&scanline, tile) {
                        self.reveal(tile);
                    }
                    if let Some(prev_tile) = prev_tile {
                        if self.is_opaque(prev_tile) && !self.is_opaque(tile) {
                            scanline.start_slope = slope(tile);
                        }
                        if !self.is_opaque(prev_tile) && self.is_opaque(tile) {
                            let mut next_line = scanline.next();
                            next_line.end_slope = slope(tile);
                            stack.push(next_line);
                        }
                    }
                    prev_tile = Some(tile);
                }
            }
            if let Some(prev_tile) = prev_tile {
                if !self.is_opaque(prev_tile) {
                    stack.push(scanline.next());
                }
            }
        }
    }
}

/// Calculates field-of-view for a map that supports Algorithm2D, returning a HashSet. This is a bit faster
/// than coercing the results into a vector, since internally it uses the set for de-duplication.
pub fn field_of_view_set<T: Algorithm2D>(
    origin: Point,
    radius: i32,
    fov_check: &T,
) -> HashSet<Point> {
    // NOTE: Symmetric shadowcasting algorithm adapted from https://www.albertford.com/shadowcasting/ (CC0)

    let mut visible_points: HashSet<Point> = HashSet::with_capacity((4 * radius * radius) as usize);
    visible_points.insert(origin);

    let radius_plus_half = Rational32::from_integer(radius) + Rational32::new(1, 2);
    let radius_plus_half_2 = radius_plus_half * radius_plus_half;
    let radius_2 = radius * radius;

    for cardinal in [
        Cardinal::North,
        Cardinal::East,
        Cardinal::South,
        Cardinal::West,
    ] {
        let mut scanner = FovScanner {
            radius_2,
            radius_plus_half_2,
            quadrant: Quadrant::new(cardinal, origin),
            fov_check,
            visible_points: &mut visible_points,
        };
        let first_line = Scanline::with_integers(1, -1, 1);
        scanner.scan(first_line);
    }

    visible_points
        .iter()
        .copied()
        .filter(|p| fov_check.in_bounds(*p))
        .collect()
}

/// Calculates field-of-view for a map that supports Algorithm2D.
pub fn field_of_view<T: Algorithm2D>(start: Point, range: i32, fov_check: &T) -> Vec<Point> {
    field_of_view_set(start, range, fov_check)
        .into_iter()
        .collect()
}

#[cfg(test)]
mod tests {

    use crate::prelude::*;
    use bracket_algorithm_traits::prelude::{Algorithm2D, BaseMap};
    use bracket_geometry::prelude::{BresenhamCircle, Point};
    use std::cmp::max;
    use std::collections::HashSet;
    use std::hash::Hash;

    const TESTMAP_W: usize = 20;
    const TESTMAP_H: usize = 20;
    const TESTMAP_TILES: usize = (TESTMAP_W * TESTMAP_H) as usize;

    struct Map {
        pub tiles: Vec<bool>,
    }

    impl Map {
        fn new() -> Map {
            Map {
                tiles: vec![false; TESTMAP_TILES],
            }
        }
    }

    // The map needs to be see-through for the tests to check FOV
    impl BaseMap for Map {
        fn is_opaque(&self, idx: usize) -> bool {
            self.tiles[idx]
        }
    }

    impl Algorithm2D for Map {
        fn dimensions(&self) -> Point {
            Point::new(TESTMAP_W, TESTMAP_H)
        }
    }

    fn has_unique_elements<T>(iter: T) -> bool
    where
        T: IntoIterator,
        T::Item: Eq + Hash,
    {
        let mut uniq = HashSet::new();
        iter.into_iter().all(move |x| uniq.insert(x))
    }

    // Tests that we are correctly de-duplicating field of view
    #[test]
    fn fov_dupes() {
        let map = Map::new();

        let visible = field_of_view(Point::new(10, 10), 8, &map);

        assert!(has_unique_elements(&visible));
    }

    // Tests that the bounds-checking trait is applying properly to field-of-view checks
    #[test]
    fn fov_bounds_check() {
        let map = Map::new();

        let visible = field_of_view(Point::new(2, 2), 8, &map);

        for t in visible.iter() {
            assert!(t.x >= 0);
            assert!(t.x < TESTMAP_W as i32 - 1);
            assert!(t.y >= 0);
            assert!(t.y < TESTMAP_H as i32 - 1);
        }
    }

    // Tests that the FOV scan does not miss any interior points
    #[test]
    fn fov_inclusive() {
        for radius in 4..=9 {
            let map = Map::new();
            let dimensions = map.dimensions();
            let c = Point::new(10, 10);
            let visible = field_of_view(c, radius, &map);
            // let max_radius_sq: i32 = visible.iter().fold(0, |max_r2, p| {
            let max_radius_sq: i32 = BresenhamCircle::new(c, radius).fold(0, |max_r2, p| {
                let r2 = (p.x - c.x) * (p.x - c.x) + (p.y - c.y) * (p.y - c.y);
                max(r2, max_r2)
            });
            /*
            for y in 0..dimensions.y {
                let mut s = "".to_string();
                for x in 0..dimensions.x {
                    let point = Point::new(x, y);
                    let c = if visible.contains(&point) {
                        '.'
                    } else {
                        '#'
                    };
                    s.push(c);
                }
                println!("{}", s);
            }
            */
            for x in 0..dimensions.x {
                for y in 0..dimensions.y {
                    let r2 = (x - c.x) * (x - c.x) + (y - c.y) * (y - c.y);
                    let point = Point::new(x, y);
                    assert!(
                        r2 >= max_radius_sq || visible.contains(&point),
                        format!("Interior point ({:?}) not in FOV({})", point, radius)
                    );
                }
            }
        }
    }

    #[test]
    fn fov_corridor() {
        let mut map = Map::new();
        let c = Point::new(10, 10);
        let radius: i32 = 5;

        for i in 0..20 {
            let idx = 9 * 20 + i;
            map.tiles[idx] = true;
            let idx = 11 * 20 + i;
            map.tiles[idx] = true;
        }

        let visible = field_of_view(c, radius, &map);
        for i in 1..radius * 2 - 2 {
            let pos = Point::new(c.x - radius + i, c.y);
            assert!(visible.contains(&pos));
            let pos = Point::new(c.x - radius + i, c.y - 1);
            assert!(visible.contains(&pos), format!("{:?} not in result", pos));
            let pos = Point::new(c.x - radius + i, c.y + 1);
            assert!(visible.contains(&pos));
        }
    }
}
