# lidar-gnss-mapping

This repository works on tightly fusing raw GNSS measuremtns, IMU with LiDAR information for ***localization and mapping***.
Currently, gnss is aligned by using IMU， which means gnss cannot be utilized in the sytem if imu was not used or went wrong.

- [x] online alignment between LiDAR and IMU
- [x] online alignment between IMU and GNSS
- [x] state estimation, localization
- [x] Local Mapping and global Mapping


***More code is coming...***  

- [ ] online alignment between LiDAR and GNSS
- [ ] reduce performance time when processing 3d point cloud
# Acknowledge
The code is adapted from [Lio](https://github.com/hyye/lio-mapping) and [GVINS](https://github.com/HKUST-Aerial-Robotics/GVINS)
