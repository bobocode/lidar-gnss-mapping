# lidar-gnss-mapping
<<<<<<< HEAD
This repository works on tightly fusing raw GNSS measuremtns with LiDAR information for localizationa and mapping. The LOAM/Lidar Odometry part is adapted and refactored from [ALOAM](https://github.com/tops666/Aloam).
- [x] online alignment
- [x] tightly-couple fusion via FGO framework for localization
- [x] 3D Map
=======
This repository works on tightly fusing raw GNSS measuremtns, IMU with LiDAR information for localizationa and mapping.
Currently, gnss is aligned by using IMU， which means gnss cannot be utilized in the sytem if imu was not used or went wrong.
***More code is coming...***  
***to do list***
- [x] online alignment between LiDAR and GNSS
- [x] reduce performance time when processing 3d point cloud
# Acknowledge
The code is adapted from [Lio](https://github.com/hyye/lio-mapping) and [GVINS](https://github.com/HKUST-Aerial-Robotics/GVINS)

>>>>>>> cc4373a4c928dc9b6acbb18a94c4c0b13b0b05a8
