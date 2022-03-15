#ifndef _ESTIMATOR_H_
#define _ESTIMATOR_H_

#include <typeinfo>
#include <stdio.h>
#include <cmath>
#include <ceres/ceres.h>
#include <opencv2/core/eigen.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_spp.hpp>
#include <std_msgs/Header.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/search/impl/flann_search.hpp>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <eigen3/Eigen/Dense>
#include <lidar_odometry_mapping/lidar_odometry_mapping.h>
#include "utility.h"
#include "parameters.h"
#include "gnss_factor.h"
#include "imu_factor.h"
#include "gnss_initializer.h"
#include "imu_initializer.h"
#include "feature.h"
#include "cloud_process.h"
#include "laser_factor.h"


//typedef pcl::PointCloud<pcl::PointXYZI>::Ptr PointCloudPtr;

class Estimator
{
public :
 Estimator();
 ~Estimator();

 void processCloud(double t, const PointCloudT &laserCloudIn);
 void processIMU(double t, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
 void processGNSS(const std::vector<ObsPtr> &gnss_mea);
 void inputEphem(EphemBasePtr ephem_ptr);
 void inputIonoParams(double ts, const std::vector<double> &iono_params);
 void inputGNSSTimeDiff(const double t_diff);

 
 void buildLocalMap(std::vector<FeaturePerFrame> &feature_frames);
 void processOdom();
 void processMapping();
 void solveOptimization();

 bool runIMUInitialization();
 
 // internal
 void clearState();

 void updateGNSSStatistics();
  
 void vector2double();
 void double2vector();
 bool failureDetection();

 Eigen::Matrix4d transform_lb;
 Eigen::Matrix3d R_WI;

 std_msgs::Header Headers[(WINDOW_SIZE + 1)];
 Eigen::Vector3d Ps[(WINDOW_SIZE + 1)];
 Eigen::Vector3d Vs[(WINDOW_SIZE + 1)];
 Eigen::Matrix3d Rs[(WINDOW_SIZE + 1)];
 Eigen::Vector3d Bas[(WINDOW_SIZE + 1)];
 Eigen::Vector3d Bgs[(WINDOW_SIZE + 1)];

 Eigen::Matrix3d pivot_r;
 Eigen::Vector3d pivot_p;

 bool first_imu, init_local_map, first_optimization;
 int frame_count;

 std::vector<FeaturePerFrame> feature_frames[(WINDOW_SIZE + 1)];

 std::vector<PointCloudPtr> corner_stack[(WINDOW_SIZE + 1)];
 std::vector<PointCloudPtr> surf_stack[(WINDOW_SIZE + 1)];
 std::vector<PointCloudPtr> corner_less_stack[(WINDOW_SIZE + 1)];
 std::vector<PointCloudPtr> surf_less_stack[(WINDOW_SIZE + 1)];

 std::shared_ptr<LidarOdometryMapping> point_process;
 std::shared_ptr<CloudProcess> cloud_process;

 std::vector<double, PairTimeLaserTransform> all_laser_transforms;
 std::vector<Eigen::Matrix4d> local_transforms;

 PointCloudPtr local_surf_points, local_surf_filtered_points;
 PointCloudPtr local_corner_points, local_corner_filtered_points;

 PointCloudPtr localMap;
 PointCloudPtr globalMap;

 // kd-tree
  pcl::KdTreeFLANN<PointType>::Ptr kdtreeCornerFromMap;
  pcl::KdTreeFLANN<PointType>::Ptr kdtreeSurfFromMap;

  pcl::VoxelGrid<PointType> down_size_filter_corner;
  pcl::VoxelGrid<PointType> down_size_filter_surf;
  pcl::VoxelGrid<PointType> down_size_filter_map;
  
  IntegrationBase *pre_integrations[(WINDOW_SIZE + 1)];
  IntegrationBase *tmp_pre_integration;
  bool imu_enable;

  Eigen::Vector3d acc_0, gyr_0;
  Eigen::Vector3d g_vec;

  std::vector<double> dt_buf[(WINDOW_SIZE + 1)];
  std::vector<Vector3d> linear_acceleration_buf[(WINDOW_SIZE + 1)];
  std::vector<Vector3d> angular_velocity_buf[(WINDOW_SIZE + 1)];


  // GNSS related
  bool GNSSLiDARAlign();

  // GNSS related
  bool gnss_ready;
  Eigen::Vector3d anc_ecef;
  Eigen::Matrix3d R_ecef_enu;
  double yaw_enu_local;
  std::vector<ObsPtr> gnss_meas_buf[(WINDOW_SIZE+1)];
  std::vector<EphemBasePtr> gnss_ephem_buf[(WINDOW_SIZE+1)];
  std::vector<double> latest_gnss_iono_params;
  std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem;
  std::map<uint32_t, std::map<double, size_t>> sat2time_index;
  std::map<uint32_t, uint32_t> sat_track_status;
  double para_anc_ecef[3];
  double para_yaw_enu_local[1];
  double para_rcv_dt[(WINDOW_SIZE+1)*4];
  double para_rcv_ddt[WINDOW_SIZE+1];
  // GNSS statistics
  double diff_t_gnss_local;
  Eigen::Matrix3d R_enu_local;
  Eigen::Vector3d ecef_pos, enu_pos, enu_vel, enu_ypr;


 private:
  double **para_pose;
  double **para_speed_bias;
  double para_ex_pose[7];

};

#endif