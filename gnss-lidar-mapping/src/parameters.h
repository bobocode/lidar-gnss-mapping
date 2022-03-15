#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <cstdlib>
#include <ros/ros.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <fstream>

const int WINDOW_SIZE = 15;
const int OPT_WINDOW_SIZE = 5;

extern int MAX_ITERATIONS;


extern double CORNER_FILTER_SIZE;
extern double SURF_FILTER_SIZE;
extern double MAP_FILTER_SIZE;

extern std::string LIDAR_TOPIC;
extern std::string CLOUD_SURF_TOPIC;
extern std::string CLOUD_FLAT_TOPIC;
extern std::string CLOUD_SEG_TOPIC;
extern std::string CLOUD_GROUND_TOPIC;

extern bool IMU_ENABLE;
extern double BIAS_ACC_THRESHOLD;
extern double BIAS_GYR_THRESHOLD;
extern double SOLVER_TIME;
extern int NUM_ITERATIONS;
extern std::string FACTOR_GRAPH_RESULT_PATH;
extern std::string IMU_TOPIC;
extern double MSG_TIME_DELAY;

extern double G_NORM;
extern Eigen::Vector3d G;
extern Eigen::Matrix4d TRANSFORM_LB;

extern bool GNSS_ENABLE;
extern std::string GNSS_EPHEM_TOPIC;
extern std::string GNSS_GLO_EPHEM_TOPIC;
extern std::string GNSS_MEAS_TOPIC;
extern std::string GNSS_IONO_PARAMS_TOPIC;
extern std::string GNSS_TP_INFO_TOPIC;
extern std::vector<double> GNSS_IONO_DEFAULT_PARAMS;
extern bool GNSS_LOCAL_ONLINE_SYNC;
extern std::string LOCAL_TRIGGER_INFO_TOPIC;
extern double GNSS_LOCAL_TIME_DIFF;
extern double GNSS_ELEVATION_THRES;
extern double GNSS_PSR_STD_THRES;
extern double GNSS_DOPP_STD_THRES;
extern uint32_t GNSS_TRACK_NUM_THRES;
extern double GNSS_DDT_WEIGHT;
extern std::string GNSS_RESULT_PATH;
extern double MIN_SQ_DIS;
extern double MIN_PLANE_DIS;
extern double MIN_EDGE_DIS;

void readParameters(ros::NodeHandle &n);

enum StateOrder
{
 O_P = 0,
 O_R = 3,
 O_V = 6,
 O_BA = 9,
 O_BG = 12
};

enum NoiseOrder
{
 O_AN = 0,
 O_GN = 3,
 O_AW = 6,
 O_GW = 9
};


#endif