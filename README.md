# lidar-gnss-mapping

This repository works on tightly fusing raw GNSS measuremtns, IMU with LiDAR information for ***localization and mapping***.
Currently, gnss is aligned by using IMU， which means gnss cannot be utilized in the sytem if imu was not used or went wrong.
***More code is coming...***  
***to do list***
- [x] online alignment between LiDAR and GNSS
- [x] reduce performance time when processing 3d point cloud
# Acknowledge
The code is adapted from [Lio](https://github.com/hyye/lio-mapping) and [GVINS](https://github.com/HKUST-Aerial-Robotics/GVINS)
