#include <gloam_tf/gloam_tf.h>

namespace gloam
{
void GloamTF::init()
{
  tf_buffer_ptr = std::make_shared<tf2_ros::Buffer>();
  tf_listener = std::make_shared<tf2_ros::TransformListener>(*tf_buffer_ptr);

  world_to_base_last.setOrigin(tf2::Vector3(0.0, 0.0, 0.0));
  world_to_base_last.setRotation(tf2::Quaternion(0.0, 0.0, 0.0, 1.0));
  world_stamped_LOtf_base.header.frame_id = "map";
  world_stamped_LOtf_base.child_frame_id = "base";
  world_stamped_LOtf_base.transform = tf2::toMsg(world_to_base_last);
}

void GloamTF::LO2VeloAndBase(const tf2::Transform &vel_curr_T_vel_last)
{
  // get T_world^curr = T_last^curr * T_world^last

  velo_last_to_velo_curr =  vel_curr_T_vel_last.inverse();
  geometry_msgs::Transform temp = tf2::toMsg(vel_curr_T_vel_last);  // TODO: check better solution

  if (!std::isnan(temp.translation.x) and !std::isnan(temp.translation.y) and !std::isnan(temp.translation.z) and
      !std::isnan(temp.rotation.x) and !std::isnan(temp.rotation.y) and !std::isnan(temp.rotation.z) and
      !std::isnan(temp.rotation.w))                  // avoid nan at the first couple steps
    world_to_base_last *= vel_curr_T_vel_last;  // after update, last becomes the curr
  world_stamped_LOtf_base.header.stamp = ros::Time::now();
  world_stamped_LOtf_base.transform = tf2::toMsg(world_to_base_last);
}

void GloamTF::LO2VeloStartFrame(FILE *LOFilePtr, const int count,const double time)
{
  if (count < 0)
    return;

  vel_init_to_vel_last = world_to_base_last;

  if (count == 0)
    vel_start_to_vel_last = vel_init_to_vel_last;

  vel_start_to_vel_last = vel_init_to_vel_start.inverse() * vel_init_to_vel_last;
  

  vel_start_eigen_to_vel_last = tf2::transformToEigen(tf2::toMsg(vel_start_to_vel_last)).cast<float>();


  Eigen::Vector3d tmp_T(vel_start_eigen_to_vel_last(0, 3),vel_start_eigen_to_vel_last(1, 3),vel_start_eigen_to_vel_last(2, 3));
  Eigen::Matrix3d tmp_R;
  tmp_R << vel_start_eigen_to_vel_last(0, 0),vel_start_eigen_to_vel_last(0, 1),vel_start_eigen_to_vel_last(0, 2),
                           vel_start_eigen_to_vel_last(1, 0),vel_start_eigen_to_vel_last(1, 1),vel_start_eigen_to_vel_last(1, 2),
                           vel_start_eigen_to_vel_last(2, 0),vel_start_eigen_to_vel_last(2, 1),vel_start_eigen_to_vel_last(2, 2);

  Eigen::Quaterniond tmp_Q(tmp_R);
  
  if (LOFilePtr != NULL)
  {
    fprintf(LOFilePtr, "%f %f %f %f %f %f %f %f\n", time, tmp_T.x(),tmp_T.y(),tmp_T.z(),tmp_Q.x(),tmp_Q.y(),tmp_Q.z(),tmp_Q.w());
  }

}

}  // namespace gloam