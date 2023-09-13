use super::*;

/// Adjacency matrix structure for clustering analysis.
pub struct AdjacencyMatrix{ am: Vec<Vec<bool>> , n: usize }

impl fmt::Display for AdjacencyMatrix {
    /// Pretty matrix representations.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        const MAX_RENDER: usize = 40; // Max number of rows / columns to render

        if self.n <= MAX_RENDER {
            let mut display_string = String::with_capacity(2*self.n*(self.n + 7) + 19);
            {   // Upper border of pretty matrix representation
                display_string.push_str("┌  ");
                for _ in 0..self.n { display_string.push_str("  "); }
                display_string.push_str(" ┐\n")
            }

            // Main body of matrix
            for i in 0..self.n {
                display_string.push_str("│  ");
                for j in 0..self.n {
                    if self.check_adjacency(i, j).unwrap() { display_string.push_str("X "); }
                    else { display_string.push_str("- "); }
                }
                display_string.push_str(" │\n");
            }
            {   // Lower border of pretty matrix representation
                display_string.push_str("└  ");
                for _ in 0..self.n { display_string.push_str("  "); }
                display_string.push_str(" ┘")
            }
            write!(f, "{display_string}")
        } else {

            // The length of such a display string is fixed with respect to to MAX_RENDER.
            // It is computed at compile-time in case the value of MAX_RENDER is altered in future.
            const STRING_LEN: usize = MAX_RENDER * (2*MAX_RENDER + 21) + 38;
            let mut display_string = String::with_capacity(STRING_LEN);

            {   // Upper border of pretty matrix representation
                display_string.push_str("┌  ");
                for _ in 0..MAX_RENDER { display_string.push_str("  "); }
                display_string.push_str("    ┐\n")
            }

            // Main body of matrix, up to skipped row(s)
            for i in 0..(MAX_RENDER-1) {
                display_string.push_str("│  ");
                for j in 0..(MAX_RENDER-1) {
                    if self.check_adjacency(i, j).unwrap() { display_string.push_str("X "); }
                    else { display_string.push_str("- "); }
                }
                if self.check_adjacency(i, self.n-1).unwrap() { display_string.push_str("·· X │\n"); }
                else { display_string.push_str("·· -  │\n"); }
            }

            {   // Vertical dots to indicate skipped row(s)
                display_string.push_str("│  ");
                for _ in 0..(MAX_RENDER-1) { display_string.push_str(": "); }
                display_string.push_str("   :  │\n│  ");

                // Bottom row of matrix
                for j in 0..(MAX_RENDER-1) {
                    if self.check_adjacency(self.n-1, j).unwrap() { display_string.push_str("X "); }
                    else { display_string.push_str("- "); }
                }
                display_string.push_str("·· X  │\n└  ");

                // Lower border of pretty matrix representation
                for _ in 0..MAX_RENDER { display_string.push_str("  "); }
                display_string.push_str("    ┘")
            }
            write!(f, "{display_string}")
        }
    }
}

impl AdjacencyMatrix {
    /// Build an adjacency matrix from a set of objects which implement the [`Coords`](super::Coords) trait.
    /// Distinct objects within `cutoff` of each other are marked as adjacent.
    /// If two supplied objects are coincident, they *will* be marked as adjacent as the code does not check for duplicates.
    pub fn build_from_coords<T: Coords>(coords: &Vec<T>, cutoff: &f64, opt_pbc: &Option<[f32; 3]>) -> Self {
        let cutoff_sq = cutoff * cutoff;
        let n = coords.len();
        let mut am = Vec::with_capacity(n - 1);

        for i in 0..(n-1) {
            let mut am_row = Vec::with_capacity(n - i - 1);
            for j in (i+1)..n {
                am_row.push(coords[i].dsq(&coords[j], opt_pbc) <= cutoff_sq);
            }
            am.push(am_row);
        }
        
        Self { am, n }
    }

    /// Given two indices, check if the corresponding objects are adjacent.
    pub fn check_adjacency(&self, i: usize, j: usize) -> Result<bool, &'static str> {
        if i > self.n || j > self.n { Err("index out of bounds of contact matrix") }
        else if i==j { Ok(false) }
        else if i > j { Ok(self.am[j][i-j-1]) }
        else { Ok(self.am[i][j-i-1]) }
    }

    /// Group objects into clusters using depth first search (DFS).
    /// Outputs a list of clusters (as lists of indices) in descending size order.
    pub fn cluster(&self) -> Vec<Vec<usize>> {
        let mut clusters = Vec::new();
        let mut visited: Vec<bool> = (0..self.n).map(|_| false).collect();
        
        // Iterate through matrix rows, skipping rows which have already been added to a cluster
        for i in 0..self.n {
            if !visited[i] {
                // Build full cluster using recursive dfs function
                let mut cluster = Vec::new();
                self.dfs(i, &mut visited, &mut cluster);
                cluster.sort_unstable();
                clusters.push(cluster);
            }
        }

        clusters.sort_unstable_by(|a, b| b.len().cmp(&a.len()));
        clusters
    }

    /// Recursive DFS to find current cluster -- max recursion depth = n
    fn dfs(&self, i: usize, visited: &mut Vec<bool>, cluster: &mut Vec<usize>) {
        visited[i] = true;
        cluster.push(i);
        for j in 1..self.n {
            if !visited[j] && self.check_adjacency(i, j).unwrap() {
                self.dfs(j, visited, cluster);
            }
        }
    }
}