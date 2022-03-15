#include laser_factor.h

static double sqrt_info_static = 100;

PointDistanceFactor::PointDistanceFactor(const Eigen::Vector3d &point,
                                         const Eigen::Vector4d &coeff,
                                         const Eigen::Matrix<double, 6, 6> info_mat) : point_{point},
                                                                                       coeff_{coeff},
                                                                                       info_mat_{info_mat} {
}

bool PointDistanceFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {

  Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
  Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

  Eigen::Vector3d tlb(parameters[1][0], parameters[1][1], parameters[1][2]);
  Eigen::Quaterniond qlb(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);


  Eigen::Quaterniond Qli = Qi * qlb.conjugate();
  Eigen::Vector3d Pli = Pi - Qli * tlb;

  Eigen::Vector3d w(coeff_.x(), coeff_.y(), coeff_.z());
  double b = coeff_.w();

  double residual = (w.transpose() * (Qli * point_ + Pli) + b);

  double sqrt_info = sqrt_info_static;
  // FIXME: 100 magic number, info_mat
  residuals[0] = sqrt_info * residual;

  if (jacobians) {
    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d rlb = qlb.toRotationMatrix();

    if (jacobians[0]) {
      Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor> > jacobian_pose_i(jacobians[0]);
      Eigen::Matrix<double, 1, 6> jaco_i;

      jaco_i.leftCols<3>() = w.transpose();
      jaco_i.rightCols<3>() =
          -w.transpose() * Ri * (Utility::skewSymmetric(rlb.transpose() * point_) - Utility::skewSymmetric(rlb.transpose() * tlb));

      // FIXME: 100 magic number, info_mat
      jacobian_pose_i.setZero();
      jacobian_pose_i.leftCols<6>() = sqrt_info * jaco_i /* * info_mat_*/;
      jacobian_pose_i.rightCols<1>().setZero();
    }

    if (jacobians[1]) {
      Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor> > jacobian_pose_ex(jacobians[1]);
      jacobian_pose_ex.setZero();

      Eigen::Matrix<double, 1, 6> jaco_ex;
      jaco_ex.leftCols<3>() = -w.transpose() * Ri * rlb.transpose();
      jaco_ex.rightCols<3>() =
          w.transpose() * Ri * (Utility::skewSymmetric(rlb.transpose() * point_) - Utility::skewSymmetric(rlb.transpose() * tlb));

      // FIXME: 100 magic number, info_mat
      jacobian_pose_ex.setZero();
      jacobian_pose_ex.leftCols<6>() = sqrt_info * jaco_ex;
      jacobian_pose_ex.rightCols<1>().setZero();
    }
  }

  return true;
}

PivotPointPlaneFactor::PivotPointPlaneFactor(const Eigen::Vector3d &point,
                                             const Eigen::Vector4d &coeff) : point_{point},
                                                                             coeff_{coeff} {
}

bool PivotPointPlaneFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
 
  Eigen::Vector3d P_pivot(parameters[0][0], parameters[0][1], parameters[0][2]);
  Eigen::Quaterniond Q_pivot(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

  Eigen::Vector3d Pi(parameters[1][0], parameters[1][1], parameters[1][2]);
  Eigen::Quaterniond Qi(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

  Eigen::Vector3d tlb(parameters[2][0], parameters[2][1], parameters[2][2]);
  Eigen::Quaterniond qlb(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

//  Eigen::Vector3d tlb = transform_lb_.pos;
//  Eigen::Quaterniond qlb = transform_lb_.rot;

  Eigen::Quaterniond Qlpivot = Q_pivot * qlb.conjugate();
  Eigen::Vector3d Plpivot = P_pivot - Qlpivot * tlb;

  Eigen::Quaterniond Qli = Qi * qlb.conjugate();
  Eigen::Vector3d Pli = Pi - Qli * tlb;

  Eigen::Quaterniond Qlpi = Qlpivot.conjugate() * Qli;
  Eigen::Vector3d Plpi = Qlpivot.conjugate() * (Pli - Plpivot);

  Eigen::Vector3d w(coeff_.x(), coeff_.y(), coeff_.z());
  double b = coeff_.w();

  double residual = (w.transpose() * (Qlpi * point_ + Plpi) + b);

  double sqrt_info = sqrt_info_static;

  residuals[0] = sqrt_info * residual;

  if (jacobians) {
    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d Rp = Q_pivot.toRotationMatrix();
    Eigen::Matrix3d rlb = qlb.toRotationMatrix();

    if (jacobians[0]) {
      Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor> > jacobian_pose_pivot(jacobians[0]);
      Eigen::Matrix<double, 1, 6> jaco_pivot;

      jaco_pivot.leftCols<3>() = -w.transpose() * rlb * Rp.transpose();
      jaco_pivot.rightCols<3>() =
          w.transpose() * rlb * (Utility::skewSymmetric(Rp.transpose() * Ri * rlb.transpose() * (point_ - tlb))
              + Utility::skewSymmetric(Rp.transpose() * (Pi - P_pivot)));

      jacobian_pose_pivot.setZero();
      jacobian_pose_pivot.leftCols<6>() = sqrt_info * jaco_pivot;
      jacobian_pose_pivot.rightCols<1>().setZero();
    }

    if (jacobians[1]) {
      Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor> > jacobian_pose_i(jacobians[1]);
      Eigen::Matrix<double, 1, 6> jaco_i;

      jaco_i.leftCols<3>() = w.transpose() * rlb * Rp.transpose();
      jaco_i.rightCols<3>() =
          w.transpose() * rlb * Rp.transpose() * Ri * (-Utility::skewSymmetric(rlb.transpose() * point_)
              + Utility::skewSymmetric(rlb.transpose() * tlb));

      jacobian_pose_i.setZero();
      jacobian_pose_i.leftCols<6>() = sqrt_info * jaco_i;
      jacobian_pose_i.rightCols<1>().setZero();
    }

    if (jacobians[2]) {
      Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor> > jacobian_pose_ex(jacobians[2]);
      jacobian_pose_ex.setZero();

      Eigen::Matrix3d I3x3;
      I3x3.setIdentity();

      Eigen::Matrix3d right_info_mat;
      right_info_mat.setIdentity();
      right_info_mat(2, 2) = 1e-6;
      right_info_mat = Qlpivot.conjugate().normalized() * right_info_mat * Qlpivot.normalized();

      Eigen::Matrix<double, 1, 6> jaco_ex;
      //  NOTE: planar extrinsic
//       jaco_ex.leftCols<3>() = w.transpose() * (I3x3 - rlb * Rp.transpose() * Ri * rlb.transpose()) * right_info_mat;
      jaco_ex.leftCols<3>() = w.transpose() * (I3x3 - rlb * Rp.transpose() * Ri * rlb.transpose());
      jaco_ex.rightCols<3>() =
          w.transpose() * rlb * (-Utility::skewSymmetric(Rp.transpose() * Ri * rlb.transpose() * (point_ - tlb))
              + Rp.transpose() * Ri * Utility::skewSymmetric(rlb.transpose() * (point_ - tlb))
              - Utility::skewSymmetric(Rp.transpose() * (Pi - P_pivot)));

      jacobian_pose_ex.setZero();
      jacobian_pose_ex.leftCols<6>() = sqrt_info * jaco_ex;
      jacobian_pose_ex.rightCols<1>().setZero();
    }
  }

  return true;
}

PlaneProjectionFactor::PlaneProjectionFactor(const Eigen::Vector4d &local_coeffi,
                                             const Eigen::Vector4d &local_coeffj,
                                             const double &score)
    : local_coeffi_{local_coeffi}, local_coeffj_{local_coeffj}, score_{score} {}

bool PlaneProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
  Eigen::Vector3d Pi{parameters[0][0], parameters[0][1], parameters[0][2]};
  Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

  Eigen::Vector3d Pj{parameters[1][0], parameters[1][1], parameters[1][2]};
  Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

  Eigen::Vector3d tlb{parameters[2][0], parameters[2][1], parameters[2][2]};
  Eigen::Quaterniond qlb(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

  Eigen::Quaterniond Qli = Qi * qlb.conjugate();
  Eigen::Vector3d Pli = Pi - Qli * tlb;
  Eigen::Quaterniond Qlj = Qj * qlb.conjugate();
  Eigen::Vector3d Plj = Pj - Qlj * tlb;

  Twist<double> transform_li{Qli, Pli};
  Twist<double> transform_lj{Qlj, Plj};
  Twist<double> transform_lb{qlb, tlb};

  Eigen::Vector3d wi = local_coeffi_.head<3>();
  double bi = local_coeffi_.w();

//  DLOG(INFO) << "wi: " << wi.transpose();

  double disi = local_coeffi_.w();
  double disj = local_coeffj_.w();
  LOG_IF(INFO, disi < 0) << "disi less than zero: " << disi;
  LOG_IF(INFO, disj < 0) << "disj less than zero: " << disj;
//  Eigen::Vector3d plane_cen_i{-disi * local_coeffi_.x(), -disi * local_coeffi_.y(), -disi * local_coeffi_.z()};
//  Eigen::Vector3d plane_cen_j{-disj * local_coeffj_.x(), -disj * local_coeffj_.y(), -disj * local_coeffj_.z()};

  Eigen::Vector4d
      coeffi_in_j = (transform_li.inverse() * transform_lj).transform().matrix().transpose() * local_coeffi_;

  if (coeffi_in_j.w() < 0) {
//    LOG_IF(INFO, coeffi_in_j.w() < 0) << "disi_in_j less than zero: " << coeffi_in_j.w();

    coeffi_in_j = (-coeffi_in_j).eval();

    LOG_IF(INFO, coeffi_in_j.w() < 0) << "disi_in_j less than zero: " << coeffi_in_j.w();
  }

  double info = score_;
  // FIXME: 100 magic number, info_mat
  Eigen::Map<Eigen::Matrix<double, 4, 1>> residual(residuals);
  residual = info * (coeffi_in_j - local_coeffj_);

  if (jacobians) {
    Eigen::Matrix3d Ri = Qi.toRotationMatrix();
    Eigen::Matrix3d Rj = Qj.toRotationMatrix();
    Eigen::Matrix3d rlb = qlb.toRotationMatrix();

    if (jacobians[0]) {
      Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor> > jacobian_pose_i(jacobians[0]);
      Eigen::Matrix<double, 4, 6> jaco_i;

      jaco_i.setZero();
      jaco_i.bottomLeftCorner<1, 3>() = -wi.transpose() * rlb * Ri.transpose();
      jaco_i.topRightCorner<3, 3>() =
          -rlb * Rj.transpose() * Ri * Utility::skewSymmetric(rlb.transpose() * wi);
      jaco_i.bottomRightCorner<1, 3>() =
          wi.transpose() * rlb * Utility::skewSymmetric(Ri.transpose() * (Pj - Pi - Rj * rlb.transpose() * tlb));

      // FIXME: 100 magic number, info_mat
      jacobian_pose_i.setZero();
      jacobian_pose_i.leftCols<6>() = info * jaco_i;
      jacobian_pose_i.rightCols<1>().setZero();
    }

    if (jacobians[1]) {
      Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor> > jacobian_pose_j(jacobians[1]);
      Eigen::Matrix<double, 4, 6> jaco_j;

      jaco_j.setZero();
      jaco_j.bottomLeftCorner<1, 3>() = wi.transpose() * rlb * Rj.transpose();
      jaco_j.topRightCorner<3, 3>() =
          rlb * Utility::skewSymmetric(Rj.transpose() * Ri * rlb.transpose() * wi);
      jaco_j.bottomRightCorner<1, 3>() =
          wi.transpose() * rlb * Ri.transpose() * Rj * Utility::skewSymmetric(rlb.transpose() * tlb);

      // FIXME: 100 magic number, info_mat
      jacobian_pose_j.setZero();
      jacobian_pose_j.leftCols<6>() = info * jaco_j;
      jacobian_pose_j.rightCols<1>().setZero();
    }

    if (jacobians[2]) {
      Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor> > jacobian_pose_ex(jacobians[2]);
      Eigen::Matrix<double, 4, 6> jaco_ex;

      jaco_ex.setZero();
      jaco_ex.bottomLeftCorner<1, 3>() =
          -wi.transpose() * rlb * Ri.transpose() * (Rj - Ri) * rlb.transpose();
      jaco_ex.topRightCorner<3, 3>() = rlb * Rj.transpose() * Ri * Utility::skewSymmetric(rlb.transpose() * wi)
          - rlb * Utility::skewSymmetric(Rj.transpose() * Ri * rlb.transpose() * wi);
      jaco_ex.bottomRightCorner<1, 3>() =
          wi.transpose() * (-rlb * Ri.transpose() * (Rj - Ri) * Utility::skewSymmetric(rlb.transpose() * tlb)
              - rlb * Utility::skewSymmetric(Ri.transpose() * (Pj - Pi - (Rj - Ri) * rlb.transpose() * tlb)));

      // FIXME: 100 magic number, info_mat
      jacobian_pose_ex.setZero();
      jacobian_pose_ex.leftCols<6>() = info * jaco_ex;
      jacobian_pose_ex.rightCols<1>().setZero();
    }
  }

  return true;

} // Evaluate

