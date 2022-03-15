#ifndef IMU_INITIALIZER_H_
#define IMU_INITIALIZER_H_

#include <glog/logging.h>
#include <nav_msgs/Odometry.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <algorithm>
#include <map>
#include <utility>
#include <mutex>
#include <ceres/ceres.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <sophus/se3.hpp>

#include "integration_base.h"
#include "parameters.h"

struct LaserTransform {
  LaserTransform() {};
  LaserTransform(double laser_time, Eigen::Matrix4d laser_transform) : time{laser_time}, transform{laser_transform} {};

  double time;
  Eigen::Matrix4d transform;
  std::shared_ptr<IntegrationBase> pre_integration;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

typedef std::pair<double, LaserTransform> PairTimeLaserTransform;

class IMUInitializer
{
 public:
  static bool Initialization(std::vector<PairTimeLaserTransform> &all_laser_transforms,
                             Eigen::Vector3d* &Vs,
                             Eigen::Vector3d* &Bas,
                             Eigen::Vector3d* &Bgs,
                             Eigen::Vector3d &g,
                             Eigen::Matrix4d &transform_lb,
                             Eigen::Matrix3d &R_WI);

  

};


#endif