/*
 *  kdtree.hpp
 *  Header file for kdtree from nanoflann
 */

#pragma once
#include <nanoflann/nanoflann.hpp>
#include <vector>
#include <memory>
using namespace std;
using namespace nanoflann;

struct PointCloud
{
    vector<vector<double>> points;
    inline size_t kdtree_get_point_count() const { return points.size(); }
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        return points[idx][dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &) const
    {
        return false;
    }
};

// Define the KD-tree type
typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    3 /* Dimension */
    >
    KDTree;

typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    2 /* Dimension */
    >
    KDTree2;

typedef KDTreeSingleIndexAdaptor<
    L2_Simple_Adaptor<double, PointCloud>,
    PointCloud,
    -1 /* Dimension */
    >
    KDTreeN;
