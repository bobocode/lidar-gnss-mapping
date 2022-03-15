#ifndef _GNSS_FACTOR_H_
#define _GNSS_FACTOR_H_

#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <Eigen/Dense>
#include <ceres/ceres.h>

#define PSR_TO_DOPP_RATIO                   5
using namespace gnss_comm;

namespace gloam
{
 class GnssPsrDoppFactor : public ceres::SizedCostFunction<2, 7, 9, 7, 9, 1, 1, 1, 3>
 {
  public: 
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
   GnssPsrDoppFactor() = delete;
   GnssPsrDoppFactor(const ObsPtr &_obs, const EphemBasePtr &_ephem, std::vector<double> &_iono_paras, 
       const double _ratio);
   virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
   bool check_gradients(const std::vector<const double*> &parameters) const;
  private:
   const ObsPtr obs;
   const EphemBasePtr ephem;
   const std::vector<double> &iono_paras;
   double ratio;
   int freq_idx;
   double freq;
   Eigen::Vector3d sv_pos;
   Eigen::Vector3d sv_vel;
   double svdt, svddt, tgd;
   double pr_uura, dp_uura;
   double relative_sqrt_info;
 };

 class DdtSmoothFactor : public ceres::SizedCostFunction<1, 1, 1>
 {
  public: 
   DdtSmoothFactor(const double weight=1) : weight_(weight) {}
   virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
  private:
   double weight_;
 };

 class DtAnchorFactor : public ceres::SizedCostFunction<1, 1>
 {
  public: 
   DtAnchorFactor(const double dt_anchor_coeff=1000) : dt_anchor_coeff_(dt_anchor_coeff) {}
   virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
  private:
   double dt_anchor_coeff_;
 };

 class DtDdtFactor : public ceres::SizedCostFunction<1, 1, 1, 1, 1>
 {
  public: 
   DtDdtFactor() = delete;
   DtDdtFactor(const double delta_t_);
   virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
  private:
   double delta_t;
   double dt_info_coeff;
 };

class PoseAnchorFactor : public ceres::SizedCostFunction<6, 7>
{
 public: 
  PoseAnchorFactor() = delete;
  PoseAnchorFactor(const std::vector<double> anchor_value);
  virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
 private:
  Eigen::Matrix<double, 7, 1> _anchor_point;
  constexpr static double sqrt_info = 120;
};

}

#endif