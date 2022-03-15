#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>

#include "parameters.h"
#include "estimator.h"

using namespace std;
using namespace gnss_comm;

#define MAX_GNSS_LIDAR_DELAY 0.05

std::unique_ptr<Estimator> estimator_ptr;


double cur_time = -1;
double last_imu_t = -1;

std::mutex m_buf;
std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;
std::condition_variable con;

double time_diff_gnss_local;
double time_diff_valid;
int skip_parameter;


queue<std::vector<ObsPtr>> gnss_meas_buf;
queue<sensor_msgs::PointCloud2ConstPtr> points_buf;
queue<sensor_msgs::ImuConstPtr> imu_buf;

void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg)
{
    EphemPtr ephem = msg2ephem(ephem_msg);
    estimator_ptr->inputEphem(ephem);
}

void gnss_glo_ephem_callback(const GnssGloEphemMsgConstPtr &glo_ephem_msg)
{
    GloEphemPtr glo_ephem = msg2glo_ephem(glo_ephem_msg);
    estimator_ptr->inputEphem(glo_ephem);
}

void gnss_iono_params_callback(const StampedFloat64ArrayConstPtr &iono_msg)
{
    double ts = iono_msg->header.stamp.toSec();
    std::vector<double> iono_params;
    std::copy(iono_msg->data.begin(), iono_msg->data.end(), std::back_inserter(iono_params));
    assert(iono_params.size() == 8);
    estimator_ptr->inputIonoParams(ts, iono_params);
}

void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg)
{
    std::vector<ObsPtr> gnss_meas = msg2meas(meas_msg);

    latest_gnss_time = time2sec(gnss_meas[0]->time);

    // cerr << "gnss ts is " << std::setprecision(20) << time2sec(gnss_meas[0]->time) << endl;
    if (!time_diff_valid)   return;

    m_buf.lock();
    gnss_meas_buf.push(std::move(gnss_meas));
    m_buf.unlock();
    con.notify_one();
}

void imu_callback(const sensor_msgs::ImuConstPtr &raw_imu_msg)
{
    if(raw_imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = raw_imu_msg->header.stamp.toSec();

    m_buf.lock();
    imu_buf.push(raw_imu_msg);
    m_buf.unlock();
    con.notify_one();

}

void laser_callback(const sensor_msgs::PointCloud2ConstPtr &point_data_msg)
{
    m_buf.lock();
    points_buf.push(point_data_msg);
    m_buf.unlock();
    con.notify_one();

}



bool getMeasurements(std::vector<ObsPtr> &gnss_msg, sensor_msgs::PointCloud2ConstPtr &cloud_msg, std::vector<sensor_msgs::ImuConstPtr> &imu_msg)
{
    if(imu_buf.empty() || points_buf.empty() || (GNSS_ENABLE &&gnss_meas_buf.empty()))
    {
        return false;
    }

    double front_laser_ts = points_buf.front()->header.stamp.toSec();

    if(!(imu_buf.back()->header.stamp.toSec() > front_laser_ts))
    {
        return false;
    }

    double front_imu_ts = imu_buf.front()->header.stamp.toSec();

    while(!points_buf.empty() && front_imu_ts  > front_laser_ts)
    {
        points_buf.pop();
        front_laser_ts = points_buf.front()->header.stamp.toSec();
    }

    if(GNSS_ENABLE)
    {
        front_laser_ts += time_diff_gnss_local;
        double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);

        while(!gnss_meas_buf.empty() && front_gnss_ts < front_laser_ts - MAX_GNSS_LIDAR_DELAY)
        {
            gnss_meas_buf.pop();

            if(gnss_meas_buf.empty())
            {
                return false;
            }

            front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
        }

        if(gnss_meas_buf.empty())
        {
            return false;
        }else if(abs(front_gnss_ts - front_laser_ts) < MAX_GNSS_LIDAR_DELAY)
        {
            gnss_msg = gnss_meas_buf.front();
            gnss_meas_buf.pop();
        }
    }

    cloud_msg = points_buf.front();
    points_buf.pop();

    while(imu_buf.front()->header.stamp.toSec() < cloud_msg->header.stamp.toSec() + MSG_TIME_DELAY)
    {
        imu_msg.emplace_back(imu_buf.front());
        imu_buf.pop();
    }

    imu_msg.emplace_back(imu_buf.front());

    if(imu_msg.empty())
    {
        ROS_WARN("no imu bettwen two clouds.");
    }

    return true;

}

void process()
{
    while(true)
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr>> measurements;
        std::vector<sensor_msgs::ImuConstPtr> imu_msg;
        sensor_msgs::PointCloud2ConstPtr cloud_msg;
        std::vector<ObsPtr> gnss_msg;

        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
        {
            return getMeasurements(imu_msg, cloud_msg, gnss_msg);
        });
        lk.unlock();

        m_estimator.lock();

        double dx = 0, dy =0, dz = 0, rx = 0, ry = 0, rz = 0;
        for(auto &imu_data: imu_msg)
        {
            double t = imu_data->header.stamp.toSec();
            double cloud_ts = cloud_msg->header.stamp.toSec() + MSG_TIME_DELAY;

            if(t <= cloud_ts)
            {
                if(cur_time < 0)
                {
                    cur_time = t;
                }

                double dt = t- cur_time;
                ROS_ASSERT(dt >= 0);
                cur_time = t;

                dx = imu_data->linear_acceleration.x;
                dy = imu_data->linear_acceleration.y;
                dz = imu_data->linear_acceleration.z;
                rx = imu_data->angular_velocity.x;
                ry = imu_data->angular_velocity.y;
                rz = imu_data->angular_velocity.z;
                estimator_ptr->processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
            }else
            {
                double dt_1 = cloud_ts - cur_time;
                double dt_2 = t - cloud_ts;
                cur_time = cloud_ts;

                ROS_ASSERT(dt_1 >= 0);
                ROS_ASSERT(dt_2 >= 0);
                ROS_ASSERT(dt_1 + dt_2 > 0);
                double w1 = dt_2 / (dt_1 + dt_2);
                double w2 = dt_1 / (dt_1 + dt_2);
                dx = w1 * dx + w2 * imu_data->linear_acceleration.x;
                dy = w1 * dy + w2 * imu_data->linear_acceleration.y;
                dz = w1 * dz + w2 * imu_data->linear_acceleration.z;
                rx = w1 * rx + w2 * imu_data->angular_velocity.x;
                ry = w1 * ry + w2 * imu_data->angular_velocity.y;
                rz = w1 * rz + w2 * imu_data->angular_velocity.z;
                estimator_ptr->processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
            }
        }

        if(GNSS_ENABLE && !gnss_msg.empty())
        {
            estimator_ptr->processGNSS(gnss_msg);
        }


        pcl::PointCloud<pcl::PointXYZI> laset_cloud_in;
        pcl::fromROSMsg(*cloud_msg, laser_cloud_in);
        estimator_ptr->processCloud(cloud_msg->header.stamp.toSec(), laser_cloud_in);
        m_estimator.unlock();

    }
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "gloam_node");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n);
    
#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    next_pulse_time_valid = false;
    time_diff_valid = false;
    latest_gnss_time = -1;
    tmp_last_laser_time = -1;
    laser_msg_counter = 0;

    if (GNSS_ENABLE)
        skip_parameter = -1;
    else
        skip_parameter = 0;

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 100, imu_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_cloud = n.subscribe(LASER_TOPIC, 10, laser_callback);
    
    
    ros::Subscriber sub_ephem, sub_glo_ephem, sub_gnss_meas, sub_gnss_iono_params;
    ros::Subscriber sub_gnss_time_pluse_info, sub_local_trigger_info;

    if (GNSS_ENABLE)
    {
        sub_ephem = n.subscribe(GNSS_EPHEM_TOPIC, 100, gnss_ephem_callback);
        sub_glo_ephem = n.subscribe(GNSS_GLO_EPHEM_TOPIC, 100, gnss_glo_ephem_callback);
        sub_gnss_meas = n.subscribe(GNSS_MEAS_TOPIC, 100, gnss_meas_callback);
        sub_gnss_iono_params = n.subscribe(GNSS_IONO_PARAMS_TOPIC, 100, gnss_iono_params_callback);

        
        time_diff_gnss_local = GNSS_LOCAL_TIME_DIFF;
        estimator_ptr->inputGNSSTimeDiff(time_diff_gnss_local);
        time_diff_valid = true;
        
    }

    std::thread measurement_process{process};
    ros::spin();

    return 0;
}





