#pragma once

#include <vector>

namespace mtetcol {

class UnionFind
{
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    // Initialize n disjoint sets (0 to n-1)
    UnionFind(int n)
        : parent(n)
        , rank(n, 0)
    {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    // Find with path compression
    int find(int x)
    {
        if (parent[x] != x) parent[x] = find(parent[x]); // Path compression
        return parent[x];
    }

    // Union by rank
    bool unite(int x, int y)
    {
        int xr = find(x);
        int yr = find(y);
        if (xr == yr) return false; // Already in the same set

        if (rank[xr] < rank[yr]) {
            parent[xr] = yr;
        } else if (rank[xr] > rank[yr]) {
            parent[yr] = xr;
        } else {
            parent[yr] = xr;
            rank[xr]++;
        }
        return true;
    }

    // Check if two elements are in the same set
    bool connected(int x, int y) { return find(x) == find(y); }
};

} // namespace mtetcol
