#pragma once

#include <geometry_msgs/TransformStamped.h>
#include <tf2/LinearMath/Transform.h>
#include <tf2_eigen/tf2_eigen.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>
#include <tf2_ros/static_transform_broadcaster.h>
#include <tf2_ros/transform_broadcaster.h>
#include <tf2_ros/transform_listener.h>
// #include <tf/transform_datatypes.h>

namespace gloam
{
class GloamTF
{
public:
  void init();
  void processStaticTransform();
  void LO2VeloAndBase(const tf2::Transform &vel_curr_T_vel_last);

  void LO2VeloStartFrame(FILE *LOFilePtr, const int count,const double time);

  tf2_ros::StaticTransformBroadcaster static_broadcaster;
  tf2_ros::TransformBroadcaster dynamic_broadcaster;

  std::shared_ptr<tf2_ros::Buffer> tf_buffer_ptr;
  std::shared_ptr<tf2_ros::TransformListener> tf_listener;

  geometry_msgs::TransformStamped world_stamped_LOtf_base;
  tf2::Transform world_to_base_last, vel_init_to_vel_last, vel_init_to_vel_start, vel_start_to_vel_last,velo_last_to_velo_curr;
  Eigen::Isometry3f vel_start_eigen_to_vel_last;
  tf2::Transform base_prev_to_base_curr, vel_curr_to_vel_prev;

};
}  // namespace gloam