#include "estimator.h"

Estimator::Estimator()
{
 para_pose = new double *[WINDOW_SIZE+1];
 para_speed_bias = new double *[WINDOW_SIZE+1];

 for(int i = 0; i < WINDOW_SIZE +1; ++i)
 {
  para_pose[i] = new double[7];
  para_speed_bias[i] = new double[9];
 }

 clearState();
}

Estimator:: ~Estimator()
{
 for(int i =0; i < WINDOW_SIZE + 1; ++i)
 {
  delete[] para_pose[i];
  delete[] para_speed_bias[i];
 }

 delete[] para_pose;
 delete[] para_speed_bias;
}

void Estimator::clearState()
{
 for(int i =0; i < WINDOW_SIZE + 1; i++)
 {
  Rs[i].setIdentity();
  Ps[i].setZero();
  Vs[i].setZero();
  Bas[i].setZero();
  Bgs[i].setZero();

  dt_buf[i].clear();
  linear_acceleration_buf[i].clear();
  angular_velocity_buf[i].clear();

  corner_stack[i].reset();
  corner_less_stack[i].reset();
  surf_stack[i].reset();
  surf_less_stack[i].reset();

  feature_frames.clear();
  
  if(pre_integrations[i] != nullptr)
  {
   delete pre_integrations[i];
  }

  pre_integrations[i] = nullptr;
 }

 down_size_filter_corner.setLeafSize(CORNER_FILTER_SIZE, CORNER_FILTER_SIZE, CORNER_FILTER_SIZE);
 down_size_filter_surf.setLeafSize(SURF_FILTER_SIZE, SURF_FILTER_SIZE, SURF_FILTER_SIZE);
 down_size_filter_map.setLeafSize(MAP_FILTER_SIZE,MAP_FILTER_SIZE,MAP_FILTER_SIZE);

 R_WI.setIdentity();
 cloud_process = std::make_shared<CloudProcess>();
 point_process = std::make_shared<LidarOdometryMapping>();
 all_laser_transforms.clear();

 cloud_process->min_edge_dis = MIN_EDGE_DIS;
 cloud_process->min_plane_dis = MIN_PLANE_DIS;
 cloud_process->min_sq_dis = MIN_SQ_DIS;
 cloud_process->num_max_iterations = MAX_ITERATIONS;

 first_imu = false;
 init_local_map = false;

}

void Estimator::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{

 if(!first_imu)
 {
  first_imu = true;
  acc_0 = linear_acceleration;
  gyr_0 = angular_velocity;
 }

 if(!pre_integrations[frame_count])
 {
  pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
 }

 if(frame_count != 0)
 {
  pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
  tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);

  dt_buf[frame_count].push_back(dt);
  linear_acceleration_buf[frame_count].push_back(linear_acceleration);
  angular_velocity_buf[frame_count].push_back(angular_velocity);

  int j = frame_count;
  Eigen::Vector3d un_acc_0 = Rs[j] * (acc_0 = Bas[j]) - G;
  Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
  Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
  Eigen::Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - G;
  Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
  Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
  Vs[j] += dt * un_acc;
 }

 acc_0 = linear_acceleration;
 gyr_0 = angular_velocity;

}

void Estimator::processCloud(double t, const pcl::PointCloud<pcl::PointXYZI>& laserCloudIn)
{
 int idx = frame_count;
 point_process->scanRegistrationIO(laserCloudIn);
 corner_stack[idx] = point_process->cornerPointsSharp;
 surf_stack[idx] = point_process->surfPointsFlat;

 corner_less_stack[idx] = point_process->cornerPointsLessSharp;
 surf_less_stack[idx] = point_process->surfPointsLessFlat;

 Eigen::Matrix4d transform_in = Eigen::Matrix4d::Identity();

 transform_in.block<3,3>(0,0) = Rs[idx] * transform_lb.block<3,3>(0,0).transpose();
 transform_in.block<3,1>(0,3) = -1 * Rs[idx] * transform_lb.block<3,3>(0,0).transpose() * transform_lb.block<3,1>(0,3) + Ps[idx];

 // point_process->laserOdometryIO();
 // transform_in.block<3,3>(0,0) = point_process->q_wodom_curr.toRotationMatrix();
 // transform_in.block<3,1>(0,3) = point_process->t_wodom_curr;

 LaserTransform laser_transform(t, transform_in);
 all_laser_transforms.push_back(std::make_pair(t, laser_transform));
}

void Estimator::processOdom()
{
 if(frame_count < WINDOW_SIZE)
 {
  return;
 }

 
 

}

void Estimator::processMapping()
{

}

void Estimator::buildLocalMap(std::vector<FeaturePerFrame> &feature_frames);
{
 feature_frames.clear();
 
 local_surf_points.reset(new PointCloudT());
 local_surf_filtered_points.reset(new PointCloudT());
 
 local_corner_points.reset(new PointCloudT());
 local_corner_filtered_points.reset(new pointCloudT());

 PointCloudT local_normal;
 std::vector<Eigen::Matrix4d> local_transforms;

 int pivot_idx = WINDOW_SIZE - OPT_WINDOW_SIZE;

 Eigen::Vector3d Ps_pivot = Ps[pivot_idx];
 Eigen::Matrix3d Rs_pivot = Rs[pivot_idx];

 Eigen::Quaterniond rot_pivot(Rs_pivot * transform_lb.block<3,3>(0,0).transpose());
 Eigen::Vector3d pos_pivot = Ps_pivot - rot_pivot * transform_lb<3,1>(0,3);

 if(!init_local_map)
 {
  PointCloudT transformed_cloud_surf, tmp_cloud_surf;
  PointCloudT transformed_cloud_corner, tmp_cloud_corner;

  for(int i =0; i <= pivot_idx; i++)
  {
   Eigen::Vector3d Ps_i = Ps[i];
   Eigen::Matrix3d Rs_i = Rs[i];

   Eigen::Quaterniond rot_li(Rs_i * transform_lb.block<3,3>(0,0).transpose());
   Eigen::Vector3d pos_li = Ps_i - rot_li * transform_lb.block<3,1>(0,3);

   Eigen::Matrix4d transform_pivot_i = Eigen::Matrix4d::Identity();
   transform_pivot_i.block<3,3>(0,0) = rot_li.toRotationMatrix();
   transform_pivot_i.block<3,1>(0,3) = pos_li;

   pcl::transformPointCloud(*(surf_stack[i]), transformed_cloud_surf, transform_pivot_i);
   tmp_cloud_surf += transformed_cloud_surf;

   pcl::transformPointCloud(*(corner_stack[i]), transformed_cloud_corner, transform_pivot_i);
   tmp_cloud_corner += transformed_cloud_corner;
 
  }

  *(surf_stack[pivot_idx]) = tmp_cloud_surf;
  *(corner_stack[pivot_idx]) = tmp_cloud_corner;

  init_local_map = true;

 }

 for(int i =0; i < WINDOW_SIZE + 1； ++i)
 {
  Eigen::Vector3d Ps_i = Ps[i];
  Eigen::Matrix3d Rs_i = Rs[i];

  Eigen::Quaterniond rot_li(Rs_i * transform_lb.block<3,3>(0,0).transpose());
  Eigen::Vector3d pos_li = Ps_i - rot_li * transform_lb.block<3,1>(0,3);

  Eigen::Matrix4d transform_pivot_i = Eigen::Matrix4d::Identity();
  transform_pivot_i.block<3,3>(0,0) = rot_li.toRotationMatrix();
  transform_pivot_i.block<3,1>(0,3) = pos_li;

  local_transforms.push_back(transform_pivot_i);

  if(i < pivot_idx)
  {
   continue;
  }

  PointCloudT transformed_cloud_surf, transformed_cloud_corner;

  if(i != WINDOW_SIZE)
  {
   if(i == pivot_idx)
   {
    *local_surf_points += *(surf_stack[i]);
    *local_corner_points += *(corner_stack[i]);
    continue;
   }

   pcl::transformPointCloud(*(surf_stack[i]), transformed_cloud_surf, transform_pivot_i);
   pcl::transformPointCloud(*(corner_stack[i]), transformed_cloud_corner, transform_pivot_i);

   for(int p_idx = 0; p_idx < transformed_cloud_surf.size(); ++p_idx)
   {
    transformed_cloud_surf[p_idx].intensity = i;
   }
   *local_surf_points += transformed_cloud_surf;

   for(int p_idx = 0; p_idx < transformed_cloud_corner.size(); ++p_idx)
   {
    transformed_cloud_corner[p_idx].intensity = i;
   }
   *local_corner_points += transformed_cloud_corner;
  }
 }

 down_size_filter_surf.setInputCloud(local_surf_points);
 down_size_filter_surf.filter(*local_surf_filtered_points);

 down_size_filter_corner.setInputCloud(local_corner_points);
 down_size_filter_corner.filter(*local_corner_filtered_points);

 pcl::KdTreeFLANN<PointType>::Ptr kdtree_surf_from_map(new pcl::KdTreeFLANN<PointType>());
 kdtree_surf_from_map->setInputCloud(local_surf_filtered_points);

 pcl::KdTreeFLANN<PointType>::Ptr kdtree_corner_from_map(new pcl::KdTreeFLANN<PointType>());
 kdtree_corner_from_map->setInputCloud(local_corner_filtered_points);

 for(int idx = 0; idx < WINDOW_SIZE +1; ++idx)
 {
  FeaturePerFrame feature_per_frame;
  std::vector<std::unique_ptr<Feature>> features;

  if(idx > pivot_idx)
  {
   if(idx != WINDOW_SIZE)
   {
    cloud_process->computeFeatures(idx,kdtree_corner_from_map, local_corner_filtered_points, corner_stack[idx],
                                kdtree_surf_from_map, local_surf_filtered_points, surf_stack[idx], local_transforms[idx], features);
   }else
   {
    cloud_process->computeOdom(idx,kdtree_corner_from_map, local_corner_filtered_points, corner_stack[idx],
                                kdtree_surf_from_map, local_surf_filtered_points, surf_stack[idx],local_transforms[idx], features);
   }

  }

  feature_per_frame.features.assign(make_move_iterator(features.begin()), make_move_iterator(features.end()));
  feature_frames.push_back(std::move(feature_per_frame));

 }
}

void Estimator::solveOptimization()
{
 if(frame_count < WINDOW_SIZE)
 {
  ROS_ERROR("frame count is not enough for optimization.");
  return;
 }

 ceres::Problem problem;
 ceres::LossFunction *loss_function;
 loss_function = new ceres::CauchyLoss(1.0);

 std::vector<FeaturePerFrame> feature_frames;

 std::vector<double *> para_ids;

 buildLocalMap(feature_frames);

 if(imu_enable)
 {
  for(int i =0; i < OPT_WINDOW_SIZE +1; ++i)
  {
   ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
   problem.AddParameterBlock(para_pose_[i], 7, local_parameterization);
   problem.AddParameterBlock(para_speed_bias_[i], 9);
   para_ids.push_back(para_pose_[i]);
   para_ids.push_back(para_speed_bias_[i]);
  }
 }

 if(gnss_ready)
 {
  problem.AddParameterBlock(para_yaw_enu_local, 1);
  Eigen::Vector2d avg_hor_vel(0.0, 0.0);
  for (uint32_t i = 0; i <= OPT_WINDOW_SIZE; ++i)
  {
   avg_hor_vel += Vs[i].head<2>().cwiseAbs();
  }
            
  avg_hor_vel /= (OPT_WINDOW_SIZE+1);
        
   if (avg_hor_vel.norm() < 0.3)
   {
    // std::cerr << "velocity excitation not enough, fix yaw angle.\n";
    problem.SetParameterBlockConstant(para_yaw_enu_local);

   }
            
   for (uint32_t i = 0; i <= OPT_WINDOW_SIZE; ++i)
   {
    if (gnss_meas_buf[i].size() < 10)
     problem.SetParameterBlockConstant(para_yaw_enu_local);
   }
        
   problem.AddParameterBlock(para_anc_ecef, 3);
   // problem.SetParameterBlockConstant(para_anc_ecef);

   for (uint32_t i = 0; i <= OPT_WINDOW_SIZE; ++i)
   {
    for (uint32_t k = 0; k < 4; ++k)
      problem.AddParameterBlock(para_rcv_dt+i*4+k, 1);
      problem.AddParameterBlock(para_rcv_ddt+i, 1);
   }

 }

 {
  ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
  problem.AddParameterBlock(para_ex_pose_, 7, local_parameterization);
  para_ids.push_back(para_ex_pose_);
  problem.SetParameterBlockConstant(para_ex_pose_);
 }

 vector2double();

 if(first_optimization)
 {
  std::vector<double> anchor_value;
  for (uint32_t k = 0; k < 7; ++k)
  {
   anchor_value.push_back(para_Pose[0][k]);
  }
            
  PoseAnchorFactor *pose_anchor_factor = new PoseAnchorFactor(anchor_value);
  ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(pose_anchor_factor, NULL, para_pose[0]);
  first_optimization = false;
 }

 std::vector<ceres::internal::ResidualBlock *> res_ids_pim;

 if(imu_enable)
 {
  for(int i = 0; i < OPT_WINDOW_SIZE; ++i)
  {
   int j = i + 1;
   int opt_i = WINDOW_SIZE - OPT_WINDOW_SIZE + i;
   int opt_j = opt_i + 1;

   if(pre_integrations[opt_j]->sum_dt > 10.0)
   {
    continue;
   }

   IMUFactor *f = new IMUFactor(pre_integrations[opt_i]);
   ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(f,NULL,
                                   para_pose_[i],para_speed_bias_[i],
                                   para_pose_[j],para_speed_bias_[j]);

   res_ids_pim.push_back(res_id);
  }
 }

 std::vector<ceres::internal::ResidualBlock *> res_ids_psr;
 std::vector<ceres::internal::ResidualBlock *> res_ids_dtdt;
 std::vector<ceres::internal::ResidualBlock *> res_ids_smooth;
 
 if(gnss_ready)
 {
  for(int i =0; i < OPT_WINDOW_SIZE + 1; ++i)
  {
   int opt_i = WINDOW_SIZE - OPT_WINDOW_SIZE + i;
   const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[opt_i];
   const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[opt_i];

   for (uint32_t j = 0; j < curr_obs.size(); ++j)
   {
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);

    int lower_idx = -1;
    const double obs_local_ts = time2sec(curr_obs[j]->time) - diff_t_gnss_local;
    if (Headers[opt_i].stamp.toSec() > obs_local_ts)
    {
     lower_idx = (opt_i==0? 0 : opt_i-1);
    } 
    else
    {
     lower_idx = (opt_i==0? 0 : opt_i-1);
    }
     lower_idx = (opt_i==OPT_WINDOW_SIZE? OPT_WINDOW_SIZE-1 : opt_i);
    const double lower_ts = Headers[lower_idx].stamp.toSec();
    const double upper_ts = Headers[lower_idx+1].stamp.toSec();

     const double ts_ratio = (upper_ts-obs_local_ts) / (upper_ts-lower_ts);
    GnssPsrDoppFactor *gnss_factor = new GnssPsrDoppFactor(curr_obs[j], curr_ephem[j], latest_gnss_iono_params, ts_ratio);
    ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(gnss_factor, NULL, para_Pose[lower_idx], 
                 para_SpeedBias[lower_idx], para_Pose[lower_idx+1], para_SpeedBias[lower_idx+1],
                 para_rcv_dt+opt_i*4+sys_idx, para_rcv_ddt+opt_i, para_yaw_enu_local, para_anc_ecef);
    res_ids_psr.push_back(res_id);
   }
  }

  for(size_k = 0; k < 4; ++k)
  {
   for (uint32_t i = 0; i < OPT_WINDOW_SIZE; ++i)
   {
    int opt_i = WINDOW_SIZE - OPT_WINDOW_SIZE + i;
    const double gnss_dt = Headers[i+1].stamp.toSec() - Headers[opt_i].stamp.toSec();
    DtDdtFactor *dt_ddt_factor = new DtDdtFactor(gnss_dt);
    ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(dt_ddt_factor, NULL, para_rcv_dt+opt_i*4+k, 
                    para_rcv_dt+(opt_i+1)*4+k, para_rcv_ddt+i, para_rcv_ddt+opt_i+1);
    res_ids_dtdt.push_back(res_id);
   }
  }

  for (int i = 0; i < OPT_WINDOW_SIZE; ++i)
  {
   int opt_i = WINDOW_SIZE - OPT_WINDOW_SIZE +i;
   DdtSmoothFactor *ddt_smooth_factor = new DdtSmoothFactor(GNSS_DDT_WEIGHT);
   ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(ddt_smooth_factor, NULL, para_rcv_ddt+opt_i, para_rcv_ddt+opt_i+1);
   res_ids_smooth.push_back(res_id);
  }
  
 }

 std::vector<ceres::internal::ResidualBlock *> res_ids_proj;
 
 {
  for(int i =0; i < OPT_WINDOW_SIZE; i++)
  {
   int opt_i = WINDOW_SIZE - OPT_WINDOW_SIZE + i;
   FeaturePerFrame &feature_per_frame = feature_frames[opt_i];

   std::vector<std::unique_ptr<Feature>> &features = feature_per_frame.features;
   for (int j = 0; j < features.size(); ++j) 
   {
    PointPlaneFeature feature_j;
    features[j]->GetFeature(&feature_j);

    const double &s = feature_j.score;
    const Eigen::Vector3d &p_eigen = feature_j.point;
    const Eigen::Vector4d &coeff_eigen = feature_j.coeffs;

    Eigen::Matrix<double, 6, 6> info_mat_in;

    if(i != 0)
    {
      PivotPointPlaneFactor *f = new PivotPointPlaneFactor(p_eigen, coeff_eigen);
      ceres::internal::ResidualBlock *res_id = problem.AddResidualBlock(f, loss_function,
                                       para_pose[0],
                                       para_pose[i],
                                       para_ex_pose);

      res_ids_proj.push_back(res_id);
    }


   }
  }
 }


  ceres::Solver::Options options;

  options.linear_solver_type = ceres::DENSE_SCHUR;
  //options.num_threads = 8;
  options.trust_region_strategy_type = ceres::DOGLEG;
  options.max_num_iterations = 10;
  options.max_solver_time_in_seconds = 0.10;

  {
   //region residual before optimization
  }

  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  {
   //region residual after optimization
  }

  double2vector();


}

void Estimator::inputEphem(EphemBasePtr ephem_ptr)
{
 double toe = time2sec(ephem_ptr->toe);
 // if a new ephemeris comes
 if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
 {
  sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);
  sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
 }
}

void Estimator::inputIonoParams(double ts, const std::vector<double> &iono_params)
{
 if (iono_params.size() != 8)    return;

 // update ionosphere parameters
 latest_gnss_iono_params.clear();
 std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(latest_gnss_iono_params));
}

void Estimator::inputGNSSTimeDiff(const double t_diff)
{
 diff_t_gnss_local = t_diff;
}

void Estimator::processGNSS(const std::vector<ObsPtr> &gnss_meas)
{
 std::vector<ObsPtr> valid_meas;
 std::vector<EphemBasePtr> valid_ephems;
 for (auto obs : gnss_meas)
 {
  // filter according to system
  uint32_t sys = satsys(obs->sat, NULL);
  if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
   continue;

   // if not got cooresponding ephemeris yet
   if (sat2ephem.count(obs->sat) == 0)
    continue;
        
   if (obs->freqs.empty())    continue;       // no valid signal measurement
   int freq_idx = -1;
   L1_freq(obs, &freq_idx);
   if (freq_idx < 0)   continue;              // no L1 observation
   
   double obs_time = time2sec(obs->time);
   std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
   double ephem_time = EPH_VALID_SECONDS;
   size_t ephem_index = -1;
   for (auto ti : time2index)
   {
    if (std::abs(ti.first - obs_time) < ephem_time)
    {
     ephem_time = std::abs(ti.first - obs_time);
     ephem_index = ti.second;
    }
   }
   if (ephem_time >= EPH_VALID_SECONDS)
   {
    cerr << "ephemeris not valid anymore\n";
    continue;
   }
   const EphemBasePtr &best_ephem = sat2ephem.at(obs->sat).at(ephem_index);

   // filter by tracking status
   LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
   if (obs->psr_std[freq_idx]  > GNSS_PSR_STD_THRES ||
     obs->dopp_std[freq_idx] > GNSS_DOPP_STD_THRES)
   {
     sat_track_status[obs->sat] = 0;
     continue;
   }
   else
   {
    if (sat_track_status.count(obs->sat) == 0)
        sat_track_status[obs->sat] = 0;
    ++ sat_track_status[obs->sat];
   }

   if (sat_track_status[obs->sat] < GNSS_TRACK_NUM_THRES)
     continue;           // not being tracked for enough epochs

   // filter by elevation angle
   if (gnss_ready)
   {
    Eigen::Vector3d sat_ecef;
    if (sys == SYS_GLO)
      sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
    else
        sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
    double azel[2] = {0, M_PI/2.0};
    sat_azel(ecef_pos, sat_ecef, azel);
    if (azel[1] < GNSS_ELEVATION_THRES*M_PI/180.0)
      continue;
   }
   valid_meas.push_back(obs);
   valid_ephems.push_back(best_ephem);
  }
    
  gnss_meas_buf[frame_count] = valid_meas;
  gnss_ephem_buf[frame_count] = valid_ephems;
}

bool Estimator::GNSSLiDARAlign()
{
 if(gnss_ready)
 {
  return true;
 }

 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
 {
  if (gnss_meas_buf[i].empty() || gnss_meas_buf[i].size() < 10)
   return false;
 }

 // check horizontal velocity excitation
 Eigen::Vector2d avg_hor_vel(0.0, 0.0);
 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  avg_hor_vel += Vs[i].head<2>().cwiseAbs();
 avg_hor_vel /= (WINDOW_SIZE+1);
 if (avg_hor_vel.norm() < 0.3)
 {
  std::cerr << "velocity excitation not enough for GNSS-VI alignment.\n";
  return false;
 }

 std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
 std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
 {
  curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
  curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
 }

 GNSSInitializer gnss_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, latest_gnss_iono_params);

 / 1. get a rough global location
 Eigen::Matrix<double, 7, 1> rough_xyzt;
 rough_xyzt.setZero();
 if (!gnss_vi_initializer.coarse_localization(rough_xyzt))
 {
  std::cerr << "Fail to obtain a coarse location.\n";
  return false;
 }

 // 2. perform yaw alignment
 std::vector<Eigen::Vector3d> local_vs;
 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  local_vs.push_back(Vs[i]);
 Eigen::Vector3d rough_anchor_ecef = rough_xyzt.head<3>();
 double aligned_yaw = 0;
 double aligned_rcv_ddt = 0;
 if (!gnss_vi_initializer.yaw_alignment(local_vs, rough_anchor_ecef, aligned_yaw, aligned_rcv_ddt))
 {
  std::cerr << "Fail to align ENU and local frames.\n";
  return false;
 }
 // std::cout << "aligned_yaw is " << aligned_yaw*180.0/M_PI << '\n';

 // 3. perform anchor refinement
 std::vector<Eigen::Vector3d> local_ps;
 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  local_ps.push_back(Ps[i]);
 Eigen::Matrix<double, 7, 1> refined_xyzt;
 refined_xyzt.setZero();
 if (!gnss_vi_initializer.anchor_refinement(local_ps, aligned_yaw, 
  aligned_rcv_ddt, rough_xyzt, refined_xyzt))
 {
  std::cerr << "Fail to refine anchor point.\n";
  return false;
 }

 //std::cout << "refined anchor point is " << std::setprecision(20) << refined_xyzt.head<3>().transpose() << '\n';
 
 // restore GNSS states
 uint32_t one_observed_sys = static_cast<uint32_t>(-1);
 for (uint32_t k = 0; k < 4; ++k)
 {
  if (rough_xyzt(k+3) != 0)
  {
   one_observed_sys = k;
   break;
  }
 }
 for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
 {
  para_rcv_ddt[i] = aligned_rcv_ddt;
  for (uint32_t k = 0; k < 4; ++k)
  {
   if (rough_xyzt(k+3) == 0)
    para_rcv_dt[i*4+k] = refined_xyzt(3+one_observed_sys) + aligned_rcv_ddt * i;
   else
    para_rcv_dt[i*4+k] = refined_xyzt(3+k) + aligned_rcv_ddt * i;
  }
 }

 anc_ecef = refined_xyzt.head<3>();
 R_ecef_enu = ecef2rotation(anc_ecef);

 yaw_enu_local = aligned_yaw;

 return true;

}

bool Estimator::runIMUInitialization()
{
 PairTimeLaserTransform laser_trans_i, laser_trans_j;
 Eigen::Vector3d sum_g;

 for(size_t i =0; i< WINDOW_SIZE; ++i)
 {
  laser_trans_j = all_laser_transforms[i+1];
  double dt = laser_trans_j.second.pre_integration->sum_dt;
  Eigen::Vector3d tmp_g = laser_trans_j.second.pre_integration->delta_v;
  sum_g += tmp_g;
 }

 Eigen::Vector3d aver_g;
 aver_g = sum_g * 1.0 / WINDOW_SIZE;

 double var = 0;

 for(size_t i = 0; i < WINDOW_SIZE; ++i)
 {
  laser_trans_j = all_laser_transforms[i+1];
  double dt = laser_trans_j.second.pre_integration->sum_dt;
  Eigen::Vector3d tmp_g = laser_trans_j.second.pre_integration->delta_v;

  var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
 }

 var = sqrt(var/WINDOW_SIZE);

 DLOG(INFO) << "IMU variation: " << var;

 if (var < 0.25) {
   ROS_INFO("IMU excitation not enough!");
   return false;
 }

 Eigen::Vector3d g_vec_in_laser;
 bool intit_result = IMUInitializer::Initialization(all_laser_transforms, Vs, Bas, Bgs, g_vec_in_laser, transform_lb , R_WI);

 for(size_t i =0; i < WINDOW_SIZE + 1; ++i)
 {
  Eigen::Matrix4d trans_li = all_laser_transforms[i].second.transform;
  Eigen::Matrix4d trans_bi = trans_li * transform_lb;

  Ps[i] = trans_bi.block<3,1>(0,3);
  Rs[i] = trans_bi.block<3,3>(0,0);
 }

 Eigen::Matrix3d R0 = R_WI.transpose();

 double yaw = Utility::R2ypr(R0 * Rs[0]).x();
 R0 = (Utility::ypr2R(Eigen::Vector3d{-yaw, 0,0}) * R0).eval();

 R_WI = R0.transpose();

 g_vec = R0 * g_vec_in_laser;

 for(int i = 0; i <= WINDOW_SIZE; i++)
 {
  pre_integrations[i]->repropagate(Bas[i], Bgs[i]);
 }

 Eigen::Matrix3d rot_diff = R0;

 for(int i =0; i <= WINDOW_SIZE; i++)
 {
  Ps[i] = (rot_diff * Ps[i]).eval();
  Rs[i] = (rot_diff * Rs[i]).eval();
  Vs[i] = (rot_diff * Vs[i]).eval();
 }

 DLOG(WARNING) << "refined gravity:  " << g_vec.transpose();

 if(!intit_result)
 {
  return false;
 }else
 {
  return true;
 }


}




