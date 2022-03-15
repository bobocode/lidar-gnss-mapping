#include "cloud_process.h"

void CloudProcess::computeFeatures(const int id,
                      const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_corner_from_map,
                      const PointCloudPtr &corner_points_sharp,
                      const PointCloudPtr &corner_filtered_points,
                      const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_surf_from_map,
                      const PointCloudPtr &surf_points_flat,
                      const PointCloudPtr &surf_filtered_points, 
                      const Eigen::Matrix4d &transform_in, std::vector<unique_ptr<Feature>> &features)
{
 PointType point_sel, point_ori, coeff1, coeff2, point_proj;

 std::vector<int> point_search_idx(5, 0);
 std::vector<float> point_search_sq_dis(5, 0);
 Eigen::Matrix<float, 5, 3> mat_A0;
 Eigen::Matrix<float, 5, 1> mat_B0;
 Eigen::Vector3f mat_X0;
 Eigen::Matrix3f mat_A1;
 Eigen::Matrix<float, 1, 3> mat_D1;
 Eigen::Matrix3f mat_V1;

 mat_A0.setZero();
 mat_B0.setConstant(-1);
 mat_X0.setZero();

 mat_A1.setZero();
 mat_D1.setZero();
 mat_V1.setZero();

 PointCloudT laser_cloud_ori;
 PointCloudT coeff_sel;
 std::vector<float> scores;

 PointCloudPtr &surf_points_ori = surf_stack[id];
 size_t surf_points_size = surf_points_ori->points.size();
 
 PointCloudPtr &corner_points_ori = corner_stack[id];
 size_t corner_points_size = corner_points_ori->points.size();

 for(int i = 0; i < surf_points_size; i++)
 {
  point_ori = surf_points_ori->points[i];
  pointAssociateToMap(point_ori, point_sel, transform_in);

  int num_neighbors = 5;
  kdtree_surf_from_map->nearestKSearch(point_sel, num_neighbors, point_search_idx, point_search_sq_dis);

  if(point_search_sq_dis[num_neighbors -1] < min_sq_dis)
  {
   for(int j =0; j < num_neighbors; j++)
   {
    mat_A0(j, 0) = surf_filtered_points->points[point_search_idx[j]].x;
    mat_A0(j, 1) = surf_filtered_points->points[point_search_idx[j]].y;
    mat_A0(j, 2) = surf_filtered_points->points[point_search_idx[j]].z;
   }

   mat_X0 = mat_A0.colPivHouseholderQr().solve(mat_B0);

   float pa = mat_X0(0, 0);
   float pb = mat_X0(1, 0);
   float pc = mat_X0(2, 0);
   float pd = 1;

   float ps = sqrt(pa * pa + pb * pb + pc * pc);
   pa /= ps;
   pb /= ps;
   pc /= ps;
   pd /= ps;

    bool planeValid = true;
    for (int j = 0; j < num_neighbors; j++) 
    {
     if (fabs(pa * surf_filtered_points->points[point_search_idx[j]].x +
      pb * surf_filtered_points->points[point_search_idx[j]].y +
      pc * surf_filtered_points->points[point_search_idx[j]].z + pd) > min_plane_dis) {
      planeValid = false;
      break;
     }
    }

    if(planeValid)
    {
     float pd2 = pa * point_sel.x + pb * point_sel.y + pc * point_sel.z + pd;
     float s = 1 - 0.9f * fabs(pd2) / sqrt(CalcPointDistance(point_sel));

     coeff1.x = s * pa;
     coeff1.y = s * pb;
     coeff1.z = s * pc;
     coeff1.intensity = s * pd;

     bool is_in_laser_fov = false;
     PointType transform_pos;
     PointType point_on_z_axis;

     point_on_z_axis.x = 0.0;
     point_on_z_axis.y = 0.0;
     point_on_z_axis.z = 10.0;
     PointAssociateToMap(point_on_z_axis, point_on_z_axis, transform_in);

     transform_pos.x = transform_to_local.pos.x();
     transform_pos.y = transform_to_local.pos.y();
     transform_pos.z = transform_to_local.pos.z();
     float squared_side1 = CalcSquaredDiff(transform_pos, point_sel);
     float squared_side2 = CalcSquaredDiff(point_on_z_axis, point_sel);

     float check1 = 100.0f + squared_side1 - squared_side2
            - 10.0f * sqrt(3.0f) * sqrt(squared_side1);

     float check2 = 100.0f + squared_side1 - squared_side2
            + 10.0f * sqrt(3.0f) * sqrt(squared_side1);

     if (check1 < 0 && check2 > 0) 
     { /// within +-60 degree
          is_in_laser_fov = true;
     }

     if (s > 0.1 && is_in_laser_fov) 
     {
      unique_ptr<PointPlaneFeature> feature = std::make_unique<PointPlaneFeature>();
      feature->score = s;
      feature->point = Eigen::Vector3d{point_ori.x, point_ori.y, point_ori.z};
      feature->coeffs = Eigen::Vector4d{coeff1.x, coeff1.y, coeff1.z, coeff1.intensity};
      features.push_back(std::move(feature));
     }
    }
  }
 }

 for(int i =0; i < corner_points_size; i++)
 {
  point_ori = corner_points_sharp->points[i];
  pointAssociateToMap(point_ori, point_sel, transform_in);
  kdtree_corner_from_map->nearestKSearch(point_sel, 5, point_search_idx, point_search_sq_dis);

  if(point_search_sq_dis[i] < min_sq_dis)
  {
   Eigen::Vector3f vc(0, 0, 0);

   for (int j = 0; j < 5; j++) 
   {
    const PointType &point_sel_tmp = corner_filtered_points->points[point_search_idx[j]];
    vc.x() += point_sel_tmp.x;
    vc.y() += point_sel_tmp.y;
    vc.z() += point_sel_tmp.z;
   }

   cv /= 5;

   Eigen::Matrix3f = mat_a;
   mat_a.setZero();
   
   for (int j = 0; j < 5; j++) 
   {
    const PointT &point_sel_tmp = corner_filtered_points->points[point_search_idx[j]];
    Eigen::Vector3f a;
    a.x() = point_sel_tmp.x - vc.x();
    a.y() = point_sel_tmp.y - vc.y();
    a.z() = point_sel_tmp.z - vc.z();

    mat_a(0, 0) += a.x() * a.x();
    mat_a(1, 0) += a.x() * a.y();
    mat_a(2, 0) += a.x() * a.z();
    mat_a(1, 1) += a.y() * a.y();
    mat_a(2, 1) += a.y() * a.z();
    mat_a(2, 2) += a.z() * a.z();
   }
   mat_A1 = mat_a / 5.0;

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> esolver(mat_A1);
   mat_D1 = esolver.eigenvalues().real();
   mat_V1 = esolver.eigenvectors().real();

   if (mat_D1(0, 2) > 3 * mat_D1(0, 1)) 
   {

    float x0 = point_sel.x;
    float y0 = point_sel.y;
    float z0 = point_sel.z;
    float x1 = vc.x() + 0.1 * mat_V1(0, 2);
    float y1 = vc.y() + 0.1 * mat_V1(1, 2);
    float z1 = vc.z() + 0.1 * mat_V1(2, 2);
    float x2 = vc.x() - 0.1 * mat_V1(0, 2);
    float y2 = vc.y() - 0.1 * mat_V1(1, 2);
    float z2 = vc.z() - 0.1 * mat_V1(2, 2);

    Eigen::Vector3f X0(x0, y0, z0);
    Eigen::Vector3f X1(x1, y1, z1);
    Eigen::Vector3f X2(x2, y2, z2);

    Eigen::Vector3f a012_vec = (X0 - X1).cross(X0 - X2);

    Eigen::Vector3f normal_to_point = ((X1 - X2).cross(a012_vec)).normalized();

    Eigen::Vector3f normal_cross_point = (X1 - X2).cross(normal_to_point);

    float a012 = a012_vec.norm();

    float l12 = (X1 - X2).norm();

    float la = normal_to_point.x();
    float lb = normal_to_point.y();
    float lc = normal_to_point.z();

    float ld2 = a012 / l12;

    point_proj = point_sel;
    point_proj.x -= la * ld2;
    point_proj.y -= lb * ld2;
    point_proj.z -= lc * ld2;

    float ld_p1 = -(normal_to_point.x() * point_proj.x + normal_to_point.y() * point_proj.y
            + normal_to_point.z() * point_proj.z);
    float ld_p2 = -(normal_cross_point.x() * point_proj.x + normal_cross_point.y() * point_proj.y
            + normal_cross_point.z() * point_proj.z);

    float s = 1 - 0.9f * fabs(ld2);

    coeff1.x = s * la;
    coeff1.y = s * lb;
    coeff1.z = s * lc;
    coeff1.intensity = s * ld_p1;

    coeff2.x = s * normal_cross_point.x();
    coeff2.y = s * normal_cross_point.y();
    coeff2.z = s * normal_cross_point.z();
    coeff2.intensity = s * ld_p2;

    bool is_in_laser_fov = false;
    PointType transform_pos;
    transform_pos.x = transform_tobe_mapped_.pos.x();
    transform_pos.y = transform_tobe_mapped_.pos.y();
    transform_pos.z = transform_tobe_mapped_.pos.z();
    float squared_side1 = CalcSquaredDiff(transform_pos, point_sel);
    float squared_side2 = CalcSquaredDiff(point_on_z_axis_, point_sel);

    float check1 = 100.0f + squared_side1 - squared_side2
            - 10.0f * sqrt(3.0f) * sqrt(squared_side1);

    float check2 = 100.0f + squared_side1 - squared_side2
            + 10.0f * sqrt(3.0f) * sqrt(squared_side1);

    if (check1 < 0 && check2 > 0) 
    { /// within +-60 degree
          is_in_laser_fov = true;
    }

    if (s > 0.1 && is_in_laser_fov) 
    {
     unique_ptr<PointPlaneFeature> feature1 = std::make_unique<PointPlaneFeature>();
     feature1->score = s * 0.5;
     feature1->point = Eigen::Vector3d{point_ori.x, point_ori.y, point_ori.z};
     feature1->coeffs = Eigen::Vector4d{coeff1.x, coeff1.y, coeff1.z, coeff1.intensity} * 0.5;
     features.push_back(std::move(feature1));

     unique_ptr<PointPlaneFeature> feature2 = std::make_unique<PointPlaneFeature>();
     feature2->score = s * 0.5;
     feature2->point = Eigen::Vector3d{point_ori.x, point_ori.y, point_ori.z};
     feature2->coeffs = Eigen::Vector4d{coeff2.x, coeff2.y, coeff2.z, coeff2.intensity} * 0.5;
     features.push_back(std::move(feature2));
    }
   }
  }
 }
}

void computeOdom(const int id, const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_corner_from_map,
                      const PointCloudPtr &corner_points_sharp,
                      const PointCloudPtr &corner_filtered_points,
                      const pcl::KdTreeFLANN<PointType>::Ptr &kdtree_surf_from_map,
                      const PointCloudPtr &surf_points_flat,
                      const PointCloudPtr &surf_filtered_points, 
                      const Eigen::Matrix4d &transform_in, std::vector<unique_ptr<Feature>> &features)
{
 bool is_degenerate = false;
 for (size_t iter_count = 0; iter_count < num_max_iterations; ++iter_count)
 {
  computeFeatures(id, kdtree_corner_from_map,corner_points_sharp, corner_filtered_points,
                      kdtree_surf_from_map, surf_points_flat, surf_filtered_points, 
                      transform_in, features);

  size_t laser_cloud_sel_size = features.size();
  Eigen::Matrix<float, Eigen::Dynamic, 6> mat_A(laser_cloud_sel_size, 6);
  Eigen::Matrix<float, 6, Eigen::Dynamic> mat_At(6, laser_cloud_sel_size);
  Eigen::Matrix<float, 6, 6> matAtA;
  Eigen::VectorXf mat_B(laser_cloud_sel_size);
  Eigen::VectorXf mat_AtB;
  Eigen::VectorXf mat_X;
  Eigen::Matrix<float, 6, 6> matP;

  PointType point_sel, point_ori, coeff;
  
  Eigen::Vector3d t_local = transform_in.block<3,1>(0,3);
  Eigen::Matrix3d r_local = transform_in.block<3,3>(0,0);
  Eigen::Matrix3d rr = r_local;

  for(int i =0; i < laser_cloud_sel_size; i++)
  {
   PointPlaneFeature feature_i;
   features[i]->GetFeature(&feature_i);
   point_ori.x = feature_i.point.x();
   point_ori.y = feature_i.point.y();
   point_ori.z = feature_i.point.z();
   coeff.x = feature_i.coeffs.x();
   coeff.y = feature_i.coeffs.y();
   coeff.z = feature_i.coeffs.z();
   coeff.intensity = feature_i.coeffs.w();

   Eigen::Vector3f p(point_ori.x, point_ori.y, point_ori.z);
   Eigen::Vector3f w(coeff.x, coeff.y, coeff.z);

   Eigen::Vector3f J_r = -w.transpose() * (local_transform.rot * Utility::SkewSymmetric(p));
   Eigen::Vector3f J_t = w.transpose();

   float d2 = w.transpose() * (local_transform.rot * p + local_transform.pos) + coeff.intensity;

   mat_A(i, 0) = J_r.x();
   mat_A(i, 1) = J_r.y();
   mat_A(i, 2) = J_r.z();
   mat_A(i, 3) = J_t.x();
   mat_A(i, 4) = J_t.y();
   mat_A(i, 5) = J_t.z();
   mat_B(i, 0) = -d2;
  }

  mat_At = mat_A.transpose();
  matAtA = mat_At * mat_A;
  mat_AtB = mat_At * mat_B;
  mat_X = matAtA.colPivHouseholderQr().solve(mat_AtB);

  if(iter_count ==0)
  {
   Eigen::Matrix<float, 1, 6> mat_E;
   Eigen::Matrix<float, 6, 6> mat_V;
   Eigen::Matrix<float, 6, 6> mat_V2;

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, 6, 6>> esolver(matAtA);
   mat_E = esolver.eigenvalues().real();
   mat_V = esolver.eigenvectors().real();

   mat_V2 = mat_V;

   is_degenerate = false;
   float eignThre[6] = {100, 100, 100, 100, 100, 100};
   for (int i = 0; i < 6; ++i) 
   {
    if (mat_E(0, i) < eignThre[i]) 
    {
     for (int j = 0; j < 6; ++j) 
     {
       mat_V2(i, j) = 0;
     }
     is_degenerate = true;
     DLOG(WARNING) << "degenerate case";
     DLOG(INFO) << mat_E;
    } else 
    {
          break;
    }
   }
   matP = mat_V2 * mat_V.inverse();

  }

  if(is_degenerate)
  {
   Eigen::Matrix<float, 6, 1> matX2(mat_X);
   mat_X = matP * matX2;
  }

  t_local.x() += mat_X(3,0);
  t_local.y() += mat_X(4,0);
  t_local.z() += mat_X(5,0);

  r_local = r_local * Utility::deltaQ(Eigen::Vector3f(mat_X(0, 0), mat_X(1, 0), mat_X(2, 0)));

  if (!isfinite(t_local.x())) t_local.x() = 0.0;
  if (!isfinite(t_local.y())) t_local.y() = 0.0;
  if (!isfinite(t_local.z())) t_local.z() = 0.0;

  float delta_r = rad2deg(Eigen::Quaterniond(rr).angularDistance(Eigen::Quaterniond(r_local)));
  float delta_t = sqrt(pow(mat_X(3, 0) * 100, 2) + pow(mat_X(4, 0) * 100, 2) + pow(mat_X(5, 0) * 100, 2));

  if (delta_r < delta_r_abort && delta_t < delta_t_abort) {
    DLOG(INFO) << "CalculateLaserOdom iter_count: " << iter_count;
    break;
  }
 }
}
 
void CloudProcess::pointAssociateToMap(PointType const* const pi, PointType* const po, Eigen::Matrix4d &transform_to_local)
{
 Eigen::Matrix3d r = transform_to_map.block<3,3>(0,0);
 Eigen::Vector3d t = transform_to_map.block<3,1>(0,3);

 Eigen::Vector3d tmp_point(pi->x, pi->y, pi->z);
 Eigen::Vector3d local_point = r * tmp_point + t;
 
 po->x = local_point.x();
 po->y = local_point.y();
 po->z = local_point.z();
 po->intensity = pi->intensity;

}
void CloudProcess::pointAssociateTobeMapped(PointType const* const pi, PointType* const po, Eigen::Matrix4d &transform_to_cur)
{
 Eigen::Matrix3d r = transform_to_cur.block<3,3>(0,0);
 Eigen::Vector3d t = transform_to_cur.block<3,1>(0,3);

 Eigen::Vector3d tmp_point(pi->x, pi->y, pi->z);
 Eigen::Vector3d cur_point = r.transpose() * (tmp_point - t);

 po->x = cur_point.x();
 po->y = cur_point.y();
 po->z = cur_point.z();
 po->intensity = pi->intensity;

}
