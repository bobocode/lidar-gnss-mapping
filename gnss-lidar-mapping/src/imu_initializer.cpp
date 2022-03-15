#include "imu_initializer.h"

Eigen::MatrixXd TangentBasis(Eigen::Vector3d &g0)
{
 Eigen::Vector3d b, c;
 Eigen::Vector3d a = g0.normalized();
 Eigen::Vector3d tmp(0, 0, 1);
 if(a == tmp)
     tmp << 1, 0, 0;
 b = (tmp - a * (a.transpose() * tmp)).normalized();
 c = a.cross(b);
 Eigen::MatrixXd bc(3, 2);
 bc.block<3, 1>(0, 0) = b;
 bc.block<3, 1>(0, 1) = c;
 return bc;
}

void solveGyroscopeBias(std::vector<PairTimeLaserTransform> &all_laser_transforms, Eigen::Vector3d* Bgs)
{
 Eigen::Matrix3d A;
 Eigen::Vector3d b;
 Eigen::Vector3d delta_bg;
 A.setZero();
 b.setZero();

 for(size_t i = 0; i < all_laser_transforms.size()-1; i++)
 {
  PairTimeLaserTransform &laser_trans_i = all_laser_transforms[i];
  PairTimeLaserTransform &laser_trans_j = all_laser_transforms[i+1];

  Eigen::MatrixXd tmp_A(3,3);
  tmp_A.setZero();
  Eigen::VectorXd tmp_b(3);
  tmp_b.setZero();

  Eigen::Matrix3d Ri = laser_trans_i.second.transform.block<3,3>(0,0);
  Eigen::Matrix3d Rj = laser_trans_j.second.transform.block<3,3>(0,0);

  Eigen::Quaterniond q_ij(Ri.conjugate() * Rj);
  tmp_A = laser_trans_j.second.pre_integration->jacobian.block<3,3>(O_R, O_BG);
  tmp_b = 2 * (laser_trans_j.second.pre_integration->delta_q.conjugate() * q_ij).vec();

  A += tmp_A.transpose() * tmp_A;
  b += tmp_A.transpose() * tmp_b;

 }

 delta_bg = A.ldlt().solve(b);
 DLOG(WARNING) << "gyroscope bias initial calibration: " << delta_bg.transpose();

 for (int i = 0; i < all_laser_transforms.size(); ++i) {
   Bgs[i] += delta_bg;
 }

 for(size_t i =0; i < all_laser_transforms.size()-1; i++)
 {
  PairTimeLaserTransform &laser_trans_j = all_laser_transforms[i+1];
  laser_trans_j.second.pre_integration->repropagate(Eigen::Vector3d::Zero(), Bgs[0]);
 }

}

bool ApproximateGravity(std::vector<PairTimeLaserTransform> &all_laser_transforms, Eigen::Vector3d &g, Eigen::Matrix4d &transform_lb)
{
 size_t num_states = 3;
 size_t window_size = all_laser_transforms.size() - 1;

 if (window_size < 5) {
   DLOG(WARNING) << ">>>>>>> window size not enough <<<<<<<";
   return false;
 }

 Eigen::Matrix3d I3x3;
 I3x3.setIdentity();

 Eigen::MatrixXd A{num_states, num_states};
 A.setZero();
 Eigen::VectorXd b{num_states};
 b.setZero();

 Eigen::Vector3d &g_approx = g;

 for(size_t i =0; i < window_size-1; ++i)
 {
  PairTimeLaserTransform &laser_trans_i = all_laser_transforms[i];
  PairTimeLaserTransform &laser_trans_j = all_laser_transforms[i+1];
  PairTimeLaserTransform &laser_trans_k = all_laser_transforms[i+2];

  double dt12 = laser_trans_j.second.pre_integration->sum_dt
  double dt23 = laser_trans_k.second.pre_integration->sum_dt;

  Eigen::Vector3d dp12 = laser_trans_j.second.pre_integration->delta_p;
  Eigen::Vector3d dp23 = laser_trans_k.second.pre_integration->delta_p;
  Eigen::Vector3d dv12 = laser_trans_j.second.pre_integration->delta_v;

  Eigen::Vector3d pl1 = laser_trans_i.second.transform.block<3,1>(0,3);
  Eigen::Vector3d pl2 = laser_trans_j.second.transform.block<3,1>(0,3);
  Eigen::Vector3d pl3 = laser_trans_k.second.transform.block<3,1>(0,3);
  Eigen::Vector3d plb = transform_lb.block<3,1>(0,3);

  Eigen::Matrix3d rl1 = laser_trans_i.second.transform.block<3,3>(0,0);
  Eigen::Matrix3d rl2 = laser_trans_j.second.transform.block<3,3>(0,0);
  Eigen::Matrix3d rl3 = laser_trans_k.second.transform.block<3,3>(0,0);
  Eigen::Matrix3d rlb = transform_lb.block<3,3>(0,0);

  Eigen::MatrixXd tmp_A(3, 3);
  tmp_A.setZero();
  Eigen::VectorXd tmp_b(3);
  tmp_b.setZero();

  tmp_A = 0.5 * I3x3 * (dt12 * dt12 * dt23 + dt23 * dt23 * dt12);
  tmp_b = (pl2 - pl1) * dt23 - (pl3 - pl2) * dt12
      + (rl2 - rl1) * plb * dt23 - (rl3 - rl2) * plb * dt12
      + rl2 * rlb * dp23 * dt12 + rl1 * rlb * dv12 * dt12 * dt23
      - rl1 * rlb * dp12 * dt23;

  A += tmp_A.transpose() * tmp_A;
  b -= tmp_A.transpose() * tmp_b;

//    A += tmp_A;
//    b -= tmp_b;

  }

  A = A * 10000.0;
  b = b * 10000.0;

//  DLOG(INFO) << "A" << endl << A;
//  DLOG(INFO) << "b" << endl << b;

  g_approx = A.ldlt().solve(b);

//  g_approx = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

  DLOG(INFO) << "g_approx: " << g_approx.transpose();

  // TODO: verify g

  double g_norm = all_laser_transforms.first().second.pre_integration->config_.g_norm;

  return fabs(g_approx.norm() - g_norm) <= 1.0;

}

void RefineGravityAccBias(std::vector<PairTimeLaserTransform> &all_laser_transforms, Eigen::Vector3d* &Vs, 
                           Eigen::Vector3d* &Bgs, Eigen::Vector3d &g_approx, Eigen::Matrix4d &transform_lb, Eigen::Matrix3d &R_WI)
{
 typedef Sophus::SO3d SO3;
 Eigen::Vector3d &g_refined = g_approx;

 size_t size_velocities = all_laser_transforms.size();
 size_t num_states = size_velocities * 3 +2;

 LOG_ASSERT(size_velocities >= 5) << ">>>>>>> window size not enough <<<<<<<";

 Eigen::Matrix3d I3x3;
 I3x3.setIdentity();

 Eigen::MatrixXd A{num_states, num_states};
 A.setZero();
 Eigen::VectorXd b{num_states};
 b.setZero();
 Eigen::VectorXd x{num_states};
 x.setZero();

 double g_norm = G_NORM;
 g_refined = g_refined.normalized() * g_norm;

 g_refined = g_refined.normalized() * g_norm;

 for (int k = 0; k < 5; ++k) {

  Eigen::MatrixXd lxly(3, 2);
  lxly = TangentBasis(g_refined);

  for (size_t i = 0; i < size_velocities - 1; ++i) 
  {
   PairTimeLaserTransform &laser_trans_i = all_laser_transforms[i];
   PairTimeLaserTransform &laser_trans_j = all_laser_transforms[i + 1];

   Eigen::MatrixXd tmp_A(6, 8);
   tmp_A.setZero();
   Eigen::VectorXd tmp_b(6);
   tmp_b.setZero();

   double dt12 = laser_trans_j.second.pre_integration->sum_dt;

   Eigen::Vector3d dp12 = laser_trans_j.second.pre_integration->delta_p;
   Eigen::Vector3d dv12 = laser_trans_j.second.pre_integration->delta_v;

   Eigen::Vector3d pl1 = laser_trans_i.second.transform.block<3,1>(0,3);
   Eigen::Vector3d pl2 = laser_trans_j.second.transform.block<3,1>(0,3);
   Eigen::Vector3d plb = transform_lb.block<3,1>(0,3);

   Eigen::Matrix3d rl1 = laser_trans_i.second.transform.block<3,3>(0,0);
   Eigen::Matrix3d rl2 = laser_trans_j.second.transform.block<3,3>(0,0);
   Eigen::Matrix3d rlb = transform_lb.block<3,3>(0,0);

   tmp_A.block<3, 3>(0, 0) = dt12 * Eigen::Matrix3d::Identity();
      pl2 - pl1 - rl1 * rlb * dp12 - (rl1 - rl2) * plb - 0.5 * r_WI * gI_n * g_norm * dt12 * dt12;
   tmp_A.block<3, 2>(0, 6) = 0.5 * Eigen::Matrix3d::Identity() * lxly * dt12 * dt12;
   tmp_b.block<3, 1>(0, 0) =
       pl2 - pl1 - rl1 * rlb * dp12 - (rl1 - rl2) * plb - 0.5 * g_refined * dt12 * dt12;

   tmp_A.block<3, 3>(3, 0) = Eigen::Matrix3d::Identity();
   tmp_A.block<3, 3>(3, 3) = -Eigen::Matrix3d::Identity();

   tmp_A.block<3, 2>(3, 6) = Eigen::Matrix3d::Identity() * lxly * dt12;
   tmp_b.block<3, 1>(3, 0) = -rl1 * rlb * dv12 - g_refined * dt12;

   Eigen::Matrix<double, 6, 6> cov_inv = Eigen::Matrix<double, 6, 6>::Zero();
   
   cov_inv.setIdentity();

   Eigen::MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
   Eigen::VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

   A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
   b.segment<6>(i * 3) += r_b.head<6>();

   A.bottomRightCorner<2, 2>() += r_A.bottomRightCorner<2, 2>();
   b.tail<2>() += r_b.tail<2>();

   A.block<6, 2>(i * 3, num_states - 2) += r_A.topRightCorner<6, 2>();
   A.block<2, 6>(num_states - 2, i * 3) += r_A.bottomLeftCorner<2, 6>();

  }

  A = A * 1000.0;
  b = b * 1000.0;
  x = A.ldlt().solve(b);

  DLOG(INFO) << "k: " << k << ", x: " << x.transpose();

  // TODO: verify if the gravity is right
  Eigen::Vector2d dg = x.segment<2>(num_states - 2);
  DLOG(INFO) << "dg: " << dg.x() << ", " << dg.y();

  g_refined = (g_refined + lxly * dg).normalized() * g_norm;
 }

 Eigen::Vector3d gI_n{0.0, 0.0, -1.0};
 Eigen::Vector3d gW_n = g_refined.normalized(); // NOTE: the Lidar's world frame
 Eigen::Vector3d gIxgW = gI_n.cross(gW_n);
 Eigen::Vector3d v_WI = gIxgW / gIxgW.norm();
 double ang_WI = atan2(gIxgW.norm(), gI_n.dot(gW_n));

 Eigen::Matrix3d r_WI(SO3::exp(ang_WI * v_WI).unit_quaternion());

 R_WI = r_WI;

 // WARNING check velocity
 for (int i = 0; i < size_velocities; ++i) 
 {
   Vs[i] = x.segment<3>(i * 3);
   DLOG(INFO) << "Vs[" << i << "]" << Vs[i].transpose();
 }

}

bool IMUInitializer::Initialization(std::vector<PairTimeLaserTransform> &all_laser_transforms,
                             Eigen::Vector3d* &Vs,
                             Eigen::Vector3d* &Bas,
                             Eigen::Vector3d* &Bgs,
                             Eigen::Vector3d &g,
                             Eigen::Matrix4d &transform_lb,
                             Eigen::Matrix3d &R_WI)
{

 solveGyroscopeBias(all_laser_transforms, Bgs);
 if (!ApproximateGravity(all_laser_transforms, g, transform_lb)) {
   return false;
 };
 
 RefineGravityAccBias(all_laser_transforms, Vs, Bgs, g, transform_lb, R_WI);

 return true;

}



