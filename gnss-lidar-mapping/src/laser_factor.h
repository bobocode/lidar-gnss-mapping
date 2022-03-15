#ifndef LASER_FACTOR_H
#define LASER_FACTOR_H

#include <ceres/ceres.h>
#include <Eigen/Eigen>

#include "utility.h"

class PivotPointPlaneFactor : public ceres::SizedCostFunction<1, 7, 7, 7> {

 public:
  PivotPointPlaneFactor(const Eigen::Vector3d &point,
                        const Eigen::Vector4d &coeff);
  virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
  void Check(double **parameters);

  Eigen::Vector3d point_;
  Eigen::Vector4d coeff_;

  // TODO: necessary?
//  static Eigen::Matrix3d sqrt_info;
  static double sum_t;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

class PlaneProjectionFactor : public ceres::SizedCostFunction<4, 7, 7, 7> {

 public:
  PlaneProjectionFactor(const Eigen::Vector4d &local_coeffi, const Eigen::Vector4d &local_coeffj, const double &score);
  virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

  void Check(double **parameters);

  Eigen::Vector4d local_coeffi_;
  Eigen::Vector4d local_coeffj_;
  double score_;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

class PlaneToPlaneFactor : public ceres::SizedCostFunction<3, 7, 7, 7> {
 public:
  PlaneToPlaneFactor(const Eigen::Vector3d &pi_local, const Eigen::Vector3d &ni_local,
                     const Eigen::Vector3d &pj_local, const Eigen::Vector3d &nj_local);
  virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

  void Check(double **parameters);

  Eigen::Vector3d pi_local_;
  Eigen::Vector3d ni_local_;
  Eigen::Vector3d pj_local_;
  Eigen::Vector3d nj_local_;
//  Eigen::Matrix3d mahalanobis_;

  PointNormalFeature pfi_;
  PointNormalFeature pfj_;

  const double eps_ = 1e-8;

  static double sum_t_;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


class PointDistanceFactor : public ceres::SizedCostFunction<1, 7, 7> {

 public:
  PointDistanceFactor(const Eigen::Vector3d &point,
                      const Eigen::Vector4d &coeff,
                      const Eigen::Matrix<double, 6, 6> info_mat);
  virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
  void Check(double **parameters);

  Eigen::Vector3d point_;
  Eigen::Vector4d coeff_;
  Twist<double> transform_lb_;
  Eigen::Matrix<double, 6, 6> info_mat_;

  // TODO: necessary?
  static Eigen::Matrix3d sqrt_info;
  static double sum_t;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};


#endif