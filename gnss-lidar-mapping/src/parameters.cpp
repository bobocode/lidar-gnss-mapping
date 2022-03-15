#include "parameters.h"

std::string LIDAR_TOPIC;
std::string CLOUD_SURF_TOPIC;
std::string CLOUD_FLAT_TOPIC;
std::string CLOUD_SEG_TOPIC;
std::string CLOUD_GROUND_TOPIC;

double CORNER_FILTER_SIZE;
double SURF_FILTER_SIZE;
double MAP_FILTER_SIZE;

int MAX_ITERATIONS;

double MIN_SQ_DIS;
double MIN_PLANE_DIS;
double MIN_EDGE_DIS;

Eigen::Matrix4d TRANSFORM_LB;

double G_NORM;

bool IMU_ENABLE;
double BIAS_ACC_THRESHOLD;
double BIAS_GYR_THRESHOLD;
double SOLVER_TIME;
int NUM_ITERATIONS;
std::string FACTOR_GRAPH_RESULT_PATH;
std::string IMU_TOPIC;
double MSG_TIME_DELAY;

bool GNSS_ENABLE;
std::string GNSS_EPHEM_TOPIC;
std::string GNSS_GLO_EPHEM_TOPIC;
std::string GNSS_MEAS_TOPIC;
std::string GNSS_IONO_PARAMS_TOPIC;
std::string GNSS_TP_INFO_TOPIC;
std::vector<double> GNSS_IONO_DEFAULT_PARAMS;
bool GNSS_LOCAL_ONLINE_SYNC;
std::string LOCAL_TRIGGER_INFO_TOPIC;
double GNSS_LOCAL_TIME_DIFF;
double GNSS_ELEVATION_THRES;
double GNSS_PSR_STD_THRES;
double GNSS_DOPP_STD_THRES;
uint32_t GNSS_TRACK_NUM_THRES;
double GNSS_DDT_WEIGHT;
std::string GNSS_RESULT_PATH;

template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
 T ans;
 if (n.getParam(name, ans))
 {
  ROS_INFO_STREAM("Loaded " << name << ": " << ans);
 }
 else
 {
  ROS_ERROR_STREAM("Failed to load " << name);
  n.shutdown();
 }
 return ans;
}

void readParameters(ros::NodeHandle &n)
{
 std::string config_file;
 config_file = readParam<std::string>(n, "config_file");
 cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
 if(!fsSettings.isOpened())
 {
     std::cerr << "ERROR: Wrong path to settings" << std::endl;
 }

 fsSettings["min_sq_dis"] >> MIN_SQ_DIS;
 fsSettings["min_plane_dis"] >> MIN_PLANE_DIS;
 fsSettings["min_edge_dis"] >> MIN_EDGE_DIS;
 
 SOLVER_TIME = fsSettings["max_solver_time"];
 NUM_ITERATIONS = fsSettings["max_num_iterations"];

 std::string tmp_output_dir;
 fsSettings["output_dir"] >> tmp_output_dir;
 assert(!tmp_output_dir.empty() && "Output directory cannot be empty.\n");
 if (tmp_output_dir[0] == '~')
     tmp_output_dir.replace(0, 1, getenv("HOME"));
 char actual_output_dir[PATH_MAX+1];
 if(!realpath(tmp_output_dir.c_str(), actual_output_dir))
     std::cerr << "ERROR: Failed to obtain the real path of " << tmp_output_dir << '\n';
 std::string OUTPUT_DIR(actual_output_dir);
 FileSystemHelper::createDirectoryIfNotExists(OUTPUT_DIR.c_str());

 FACTOR_GRAPH_RESULT_PATH = OUTPUT_DIR + "/factor_graph_result.txt";
 std::ofstream fout2(FACTOR_GRAPH_RESULT_PATH, std::ios::out);
 fout2.close();

 int imu_enable_value = fsSettings["imu_enable"];
 IMU_ENABLE = (imu_enable_value == 0 ? false : true);


 ESTIMATE_EXTRINSIC = fsSettings["estimate_extrinsic"];
 if (ESTIMATE_EXTRINSIC == 2)
 {
     ROS_WARN("have no prior about extrinsic param, calibrate extrinsic param");
     RIC.push_back(Eigen::Matrix3d::Identity());
     TIC.push_back(Eigen::Vector3d::Zero());
     EX_CALIB_RESULT_PATH = OUTPUT_DIR + "/extrinsic_parameter.csv";

 }
 else 
 {
  if ( ESTIMATE_EXTRINSIC == 1)
  {
   ROS_WARN(" Optimize extrinsic param around initial guess!");
   EX_CALIB_RESULT_PATH = OUTPUT_DIR + "/extrinsic_parameter.csv";
  }
  if (ESTIMATE_EXTRINSIC == 0)
   ROS_WARN(" fix extrinsic param ");

  cv::Mat cv_R, cv_T;
  fsSettings["extrinsicRotation"] >> cv_R;
  fsSettings["extrinsicTranslation"] >> cv_T;
  Eigen::Matrix3d eigen_R;
  Eigen::Vector3d eigen_T;
  cv::cv2eigen(cv_R, eigen_R);
  cv::cv2eigen(cv_T, eigen_T);
  Eigen::Quaterniond Q(eigen_R);
  eigen_R = Q.normalized();
  RIC.push_back(eigen_R);
  TIC.push_back(eigen_T);
  ROS_INFO_STREAM("Extrinsic_R : " << std::endl << RIC[0]);
  ROS_INFO_STREAM("Extrinsic_T : " << std::endl << TIC[0].transpose());
     
 }

 int gnss_enable_value = fsSettings["gnss_enable"];
 GNSS_ENABLE = (gnss_enable_value == 0 ? false : true);

 if (GNSS_ENABLE)
 {
  fsSettings["gnss_ephem_topic"] >> GNSS_EPHEM_TOPIC;
  fsSettings["gnss_glo_ephem_topic"] >> GNSS_GLO_EPHEM_TOPIC;
  fsSettings["gnss_meas_topic"] >> GNSS_MEAS_TOPIC;
  fsSettings["gnss_iono_params_topic"] >> GNSS_IONO_PARAMS_TOPIC;
  cv::Mat cv_iono;
  fsSettings["gnss_iono_default_parameters"] >> cv_iono;
  Eigen::Matrix<double, 1, 8> eigen_iono;
  cv::cv2eigen(cv_iono, eigen_iono);
  for (uint32_t i = 0; i < 8; ++i)
   GNSS_IONO_DEFAULT_PARAMS.push_back(eigen_iono(0, i));
  
  fsSettings["gnss_tp_info_topic"] >> GNSS_TP_INFO_TOPIC;
  int gnss_local_online_sync_value = fsSettings["gnss_local_online_sync"];
  GNSS_LOCAL_ONLINE_SYNC = (gnss_local_online_sync_value == 0 ? false : true);
  if (GNSS_LOCAL_ONLINE_SYNC)
    fsSettings["local_trigger_info_topic"] >> LOCAL_TRIGGER_INFO_TOPIC;
  else
    GNSS_LOCAL_TIME_DIFF = fsSettings["gnss_local_time_diff"];

  GNSS_ELEVATION_THRES = fsSettings["gnss_elevation_thres"];
  const double gnss_ddt_sigma = fsSettings["gnss_ddt_sigma"];
  GNSS_PSR_STD_THRES = fsSettings["gnss_psr_std_thres"];
  GNSS_DOPP_STD_THRES = fsSettings["gnss_dopp_std_thres"];
  const double track_thres = fsSettings["gnss_track_num_thres"];
  GNSS_TRACK_NUM_THRES = static_cast<uint32_t>(track_thres);
  GNSS_DDT_WEIGHT = 1.0 / gnss_ddt_sigma;
  GNSS_RESULT_PATH = OUTPUT_DIR + "/gnss_result.csv";
  // clear output file
  std::ofstream gnss_output(GNSS_RESULT_PATH, std::ios::out);
  gnss_output.close();
  ROS_INFO_STREAM("GNSS enabled");
 }

 fsSettings.release();
}