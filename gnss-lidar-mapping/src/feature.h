#ifndef FEATURE_H_
#define FEATURE_H_

#include <ros/ros.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <Eigen/Eigen>
#include <ceres/ceres.h>
#include <queue>
#include <string>
#include <vector>

using std::unique_ptr;

struct Feature {
  std::string feature_name;
  virtual void GetFeature(Feature *feature) {
    DLOG(WARNING) << ">>>>>>> GetFeature not implemented <<<<<<<";
  }
};

struct PointNormalFeature : public Feature {
 public:
  PointNormalFeature() {
    feature_name = "PointNormalFeature";
  }
  PointNormalFeature(const Eigen::Vector3d &point3d_in, const Eigen::Vector3d &normal3d_in) {
    feature_name = "PointNormalFeature";
    point3d = point3d_in;
    normal3d = normal3d_in;
    diag_covariance = Eigen::Vector3d{gicp_epsilon, 1.0, 1.0}.asDiagonal();
    UpdateCovariance(normal3d);
  }

  void UpdateCovariance(const Eigen::Vector3d &normal3d_in);

  PointNormalFeature &GetFeatureInner() {
    return *this;
  }

  void GetFeature(Feature *feature) {
    PointNormalFeature *derived_feature = static_cast<PointNormalFeature *>(feature);
    *derived_feature = this->GetFeatureInner();
  }

  double gicp_epsilon = 0.001;
  Eigen::Vector3d e1{1.0, 0.0, 0.0};
  Eigen::Vector3d point3d; /// the point in the current frame
  Eigen::Vector3d normal3d; /// the normal vector in the current frame
  Eigen::Matrix3d diag_covariance;
  Eigen::Matrix3d covariance; /// the covariance in the current frame from the normal vector

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

struct PointPlaneFeature : public Feature {
  PointPlaneFeature() {
    feature_name = "PointPlaneFeature";
  }

  PointPlaneFeature(const Eigen::Vector3d &point_in, const Eigen::Vector4d &coeffs_in) {
    feature_name = "PointPlaneFeature";
    point = point_in;
    coeffs = coeffs_in;
  }

  PointPlaneFeature &GetFeatureInner() {
    return *this;
  }

  void GetFeature(Feature *feature) {
    PointPlaneFeature *derived_feature = static_cast<PointPlaneFeature *>(feature);
    *derived_feature = this->GetFeatureInner();
  }

  double score;
  Eigen::Vector3d point;
  Eigen::Vector4d coeffs;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct FeaturePerFrame {
  int id;
  std::vector<unique_ptr<Feature>> features;
};


#endif