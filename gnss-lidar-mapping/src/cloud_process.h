#ifndef CLOUD_PROCESS_H_
#define CLOUD_PROCESS_H_

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/search/impl/flann_search.hpp>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "feature.h"
#include "utility.h"

typedef pcl::PointXYZI PointType;
typedef pcl::PointCloud<PointType>::Ptr PointCloudPtr;
typedef pcl::PointCloud<PointType> PointCloudT;


class CloudProcess
{
 public:
  void computeFeatures(const int id,
                      const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_corner_from_map,
                      const PointCloudPtr &corner_points_sharp,
                      const PointCloudPtr &corner_filtered_points,
                      const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_surf_from_map,
                      const PointCloudPtr &surf_points_flat,
                      const PointCloudPtr &surf_filtered_points, 
                      const Eigen::Matrix4d &transform_in, std::vector<unique_ptr<Feature>> &features);

  void computeOdom(const int id,
                const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_corner_from_map,
                const PointCloudPtr &corner_points_sharp,
                const PointCloudPtr &corner_filtered_points,
                const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_surf_from_map,
                const PointCloudPtr &surf_points_flat,
                const PointCloudPtr &surf_filtered_points, 
                const Eigen::Matrix4d &transform_in, std::vector<unique_ptr<Feature>> &features);
 
  void pointAssociateToMap(PointType const* const pi, PointType* const po, Eigen::Matrix4d &transform_to_local);
  void pointAssociateTobeMapped(PointType const* const pi, PointType* const po, Eigen::Matrix4d &transform_to_cur);
  inline float calPointDistance(const PointType &p);

  template<typename PointType>
  inline float CalcSquaredDiff(const PointType &a, const PointType &b) {
  float diff_x = a.x - b.x;
  float diff_y = a.y - b.y;
  float diff_z = a.z - b.z;

  return diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
  }

  inline float rad2deg(float radians)
  {
    return radians * 180.0 / M_PI;
  }

  inline float deg2rad(float degrees)
  {
    return degrees * M_PI / 180.0;
  }

  template<typename PointType>
  inline float CalcSquaredDiff(const PointType &a, const PointType &b, const float &wb) {
  float diff_x = a.x - b.x * wb;
  float diff_y = a.y - b.y * wb;
  float diff_z = a.z - b.z * wb;

  return diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
  }

  template<typename PointType>
  inline float CalcPointDistance(const PointType &p) {
  return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
 }

 template<typename PointType>
  inline float CalcSquaredPointDistance(const PointType &p) {
  return p.x * p.x + p.y * p.y + p.z * p.z;
 }

public:
 double min_sq_dis;
 double min_plane_dis;
 double min_edge_dis;
 int num_max_iterations;
 double delta_r_abort;
 double delta_t_abort;
   
};



#endif